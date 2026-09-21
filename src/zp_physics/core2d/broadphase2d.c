#include "zp_physics/core2d/broadphase2d.h"
#include "zp_physics/core2d/world2d.h"
#include <string.h>
#include <assert.h>


zp_inline zp_aabb2d combine(const zp_aabb2d a, const zp_aabb2d b) {
 zp_aabb2d out;
 out._min.x = zp_min(a._min.x, b._min.x);
 out._min.y = zp_min(a._min.y, b._min.y);

 out._max.x = zp_max(a._max.x, b._max.x);
 out._max.y = zp_max(a._max.y, b._max.y);
 return out;
}


zp_inline int overlaps(zp_aabb2d a, zp_aabb2d b) {
 return a._min.x <= b._max.x &&
   a._max.x >= b._min.x &&
   a._min.y <= b._max.y &&
   a._max.y >= b._min.y;
}


zp_inline float perimeter(zp_aabb2d a) {
 float w = a._max.x - a._min.x;
 float h = a._max.y - a._min.y;
 float wh = w + h;
 return wh + wh;
}


zp_inline int contains(zp_aabb2d a, zp_aabb2d b) {
 return a._min.x <= b._min.x &&
  a._min.y <= b._min.y &&
  a._max.x >= b._max.x &&
  a._max.y >= b._max.y;
}

zp_inline zp_aabb2d expanded(zp_aabb2d a, float amount) {
 zp_aabb2d out;
 out._min.x = a._min.x - amount;
 out._min.y = a._min.y - amount;
 out._max.x = a._max.x + amount;
 out._max.y = a._max.y + amount;
 return out;
}

zp_inline int is_leaf(zp_broadphase2d_node n) {
 return n._left == ZP_POOL_NULL_ID && n._right == ZP_POOL_NULL_ID;
}


zp_noinline static void refit_upwards(zp_broadphase2d *const bp, zp_pool_id n) {
 assert(n != ZP_POOL_NULL_ID);
 zp_pool_id node = n;
 while(node != ZP_POOL_NULL_ID) {  
  zp_broadphase2d_node *node_ptr = (zp_broadphase2d_node*)zp_pool_get(&bp->_node_allocator, node);
  if(!is_leaf(*node_ptr)) {
   zp_broadphase2d_node *left = (zp_broadphase2d_node*)zp_pool_get(&bp->_node_allocator, node_ptr->_left);
   zp_broadphase2d_node *right = (zp_broadphase2d_node*)zp_pool_get(&bp->_node_allocator, node_ptr->_right);
   node_ptr->_aabb = combine(left->_aabb, right->_aabb);
  }
  node = node_ptr->_parent;
 }
}

zp_noinline static zp_pool_id find_best_sibling(zp_broadphase2d *const bp, zp_pool_id leaf_id) {
 assert(bp->_root != ZP_POOL_NULL_ID);
 
 if(is_leaf(*((zp_broadphase2d_node*)zp_pool_get(&bp->_node_allocator, bp->_root))))
  return bp->_root;
 
 zp_broadphase2d_node *leaf = (zp_broadphase2d_node*)zp_pool_get(&bp->_node_allocator, leaf_id);
 zp_broadphase2d_node *current_data;
 zp_pool_id current = bp->_root;
 do {
  zp_broadphase2d_node *node = (zp_broadphase2d_node*)zp_pool_get(&bp->_node_allocator, current);

  zp_broadphase2d_node *left = (zp_broadphase2d_node*)zp_pool_get(&bp->_node_allocator, node->_left);
  zp_broadphase2d_node *right = (zp_broadphase2d_node*)zp_pool_get(&bp->_node_allocator, node->_right);

  zp_aabb2d combined = combine(node->_aabb, leaf->_aabb);

  float inheritance_cost = 2.0f * (perimeter(combined) - perimeter(node->_aabb));

  float cost1, cost2;
  if(is_leaf(*left)) {
   zp_aabb2d box = combine(leaf->_aabb, left->_aabb);
   cost1 = perimeter(box) + inheritance_cost;
  } else {
   zp_aabb2d box = combine(leaf->_aabb, left->_aabb);

   float old_area = perimeter(left->_aabb);
   float new_area = perimeter(box);
   cost1 = (new_area - old_area) + inheritance_cost;
  }
  if(is_leaf(*right)) {
   zp_aabb2d box = combine(leaf->_aabb, right->_aabb);
   cost2 = perimeter(box) + inheritance_cost;
  } else {
   zp_aabb2d box = combine(leaf->_aabb, right->_aabb);

   float old_area = perimeter(right->_aabb);
   float new_area = perimeter(box);
   cost2 = (new_area - old_area) + inheritance_cost;
  }

  float cost = 2.0f * perimeter(combined);
  if(cost < cost1 && cost < cost2)
   break;
   
   current_data = (cost1 < cost2) ? left : right;
   current = current_data->_id;
  } while(!is_leaf(*current_data));
 return current;
}


zp_noinline static zp_pool_id create_leaf(zp_broadphase2d *const bp, const zp_aabb2d aabb, const zp_container_id body_id) {
 zp_pool_id leaf_id = zp_pool_acquire(&bp->_node_allocator);
 zp_broadphase2d_node *leaf = (zp_broadphase2d_node*)zp_pool_get(&bp->_node_allocator, leaf_id);
 leaf->_left = ZP_POOL_NULL_ID;
 leaf->_right = ZP_POOL_NULL_ID;
 leaf->_parent = ZP_POOL_NULL_ID;
 leaf->_id = leaf_id;
 leaf->_aabb = aabb;
 leaf->_body = body_id;
 return leaf_id;
}


zp_noinline static void insert_leaf(zp_broadphase2d *const bp, const zp_pool_id leaf_id) {
 if(bp->_root == ZP_POOL_NULL_ID) {
  bp->_root = leaf_id;
  return;
 }


 zp_pool_id new_parent = ZP_POOL_NULL_ID;
 zp_broadphase2d_node *new_parent_ptr = NULL;
 {
  new_parent = zp_pool_acquire(&bp->_node_allocator);
  new_parent_ptr = (zp_broadphase2d_node*)zp_pool_get(&bp->_node_allocator, new_parent);
  new_parent_ptr->_left = ZP_POOL_NULL_ID;
  new_parent_ptr->_right = ZP_POOL_NULL_ID;
  new_parent_ptr->_parent = ZP_POOL_NULL_ID;
  new_parent_ptr->_id = new_parent;
 }


 zp_broadphase2d_node *leaf = (zp_broadphase2d_node*)zp_pool_get(&bp->_node_allocator, leaf_id);
 zp_pool_id sibling = find_best_sibling(bp, leaf_id);
 zp_broadphase2d_node *sibling_ptr = (zp_broadphase2d_node*)zp_pool_get(&bp->_node_allocator, sibling);
 
 zp_pool_id old_parent = sibling_ptr->_parent;
 zp_broadphase2d_node *old_parent_ptr = NULL;
 if(old_parent != ZP_POOL_NULL_ID) 
  old_parent_ptr = (zp_broadphase2d_node*)zp_pool_get(&bp->_node_allocator, old_parent);


 new_parent_ptr->_parent = old_parent;

 new_parent_ptr->_left = sibling;
 new_parent_ptr->_right = leaf_id;

 sibling_ptr->_parent = new_parent;
 leaf->_parent = new_parent;

 if(old_parent == ZP_POOL_NULL_ID) {
  bp->_root = new_parent;
 } else {
  if(old_parent_ptr->_left == sibling)
   old_parent_ptr->_left = new_parent;
  else
   old_parent_ptr->_right = new_parent;
 }
 refit_upwards(bp, new_parent);
}

/*
 disconnect only the leaf_id, so it can be reinsert again
*/
zp_noinline static void remove_leaf(zp_broadphase2d *const bp, zp_pool_id leaf_id) {
 if(leaf_id == bp->_root) {
  bp->_root = ZP_POOL_NULL_ID;
  return;
 }
 zp_broadphase2d_node *leaf = (zp_broadphase2d_node*)zp_pool_get(&bp->_node_allocator, leaf_id);
 zp_pool_id parent = leaf->_parent;
 zp_broadphase2d_node *parent_ptr = (zp_broadphase2d_node*)zp_pool_get(&bp->_node_allocator, parent);
 
 assert(parent_ptr->_left == leaf_id || parent_ptr->_right == leaf_id);

 zp_pool_id grand_parent = parent_ptr->_parent;
 zp_broadphase2d_node *grand_parent_ptr = NULL;
 if(parent_ptr->_parent != ZP_POOL_NULL_ID)
  grand_parent_ptr = (zp_broadphase2d_node*)zp_pool_get(&bp->_node_allocator, grand_parent);

 zp_pool_id sibling = (parent_ptr->_left == leaf_id) ? parent_ptr->_right : parent_ptr->_left;

 if(parent_ptr->_parent == ZP_POOL_NULL_ID) {
  bp->_root = sibling;
  zp_broadphase2d_node *sibling_ptr = (zp_broadphase2d_node*)zp_pool_get(&bp->_node_allocator, sibling);
  sibling_ptr->_parent = ZP_POOL_NULL_ID;
 } else {
  zp_broadphase2d_node *sibling_ptr = (zp_broadphase2d_node*)zp_pool_get(&bp->_node_allocator, sibling);

  if(grand_parent_ptr->_left == parent)
   grand_parent_ptr->_left = sibling;
  else
   grand_parent_ptr->_right = sibling;
   sibling_ptr->_parent = grand_parent;
   refit_upwards(bp, grand_parent);
  }
 zp_pool_release(&bp->_node_allocator, parent);
}


zp_noinline static void traverse_pair_recursive(zp_broadphase2d *const bp, zp_pool_id a_id, zp_pool_id b_id, void (*func)(zp_container_id, zp_container_id, void *), void *ptr, uint32_t depth, zp_broadphase2d_treetranversalinfo *info) {
 if(a_id == ZP_POOL_NULL_ID || b_id == ZP_POOL_NULL_ID)
  return;

 zp_broadphase2d_node *a = (zp_broadphase2d_node *)zp_pool_get(&bp->_node_allocator, a_id);
 zp_broadphase2d_node *b = (zp_broadphase2d_node *)zp_pool_get(&bp->_node_allocator, b_id);
 
 if(!overlaps(a->_aabb, b->_aabb))
  return;

 int a_leaf = is_leaf(*a);
 int b_leaf = is_leaf(*b);
 
 info->_node_visits++;
 info->_max_depth = info->_max_depth < depth ? depth : info->_max_depth;

 if(a_leaf && b_leaf) {
  if(a_id != b_id) {
   info->_pair_test++;
   func(a->_body, b->_body, ptr);
  }
  return;
 }
 if(a_leaf) {
  traverse_pair_recursive(bp, a_id, b->_left, func, ptr, depth + 1, info);
  traverse_pair_recursive(bp, a_id, b->_right, func, ptr, depth + 1, info);
  return;
 }
 if(b_leaf) {
  traverse_pair_recursive(bp, a->_left, b_id, func, ptr, depth + 1, info);
  traverse_pair_recursive(bp, a->_right, b_id, func, ptr, depth + 1, info);
  return;
 }
 traverse_pair_recursive(bp, a->_left, b->_left, func, ptr, depth + 1, info);
 traverse_pair_recursive(bp, a->_left, b->_right, func, ptr, depth + 1, info);
 traverse_pair_recursive(bp, a->_right, b->_left, func, ptr, depth + 1, info);
 traverse_pair_recursive(bp, a->_right, b->_right, func, ptr, depth + 1, info);
}


zp_noinline static void traverse_self_recursive(zp_broadphase2d *const bp, zp_pool_id n_id, void (*func)(zp_container_id, zp_container_id, void *), void *ptr, uint32_t depth, zp_broadphase2d_treetranversalinfo *info) {
 if(n_id == ZP_POOL_NULL_ID)
  return;

 zp_broadphase2d_node *n = (zp_broadphase2d_node *)zp_pool_get(&bp->_node_allocator, n_id);
 if(is_leaf(*n)) return;

 info->_node_visits++;
 info->_max_depth = info->_max_depth < depth ? depth : info->_max_depth;

 traverse_self_recursive(bp, n->_left, func, ptr, depth + 1, info);
 traverse_self_recursive(bp, n->_right, func, ptr, depth + 1, info);
 traverse_pair_recursive(bp, n->_left, n->_right, func, ptr, depth + 1, info);
}



zp_noinline static void traverse_delete_non_leaf(zp_broadphase2d *const bp, zp_pool_id n) {
 if(n == ZP_POOL_NULL_ID) {
   return;
 }
 zp_broadphase2d_node *node = (zp_broadphase2d_node *)zp_pool_get(&bp->_node_allocator, n);

 if(is_leaf(*node))
  return;
 
 traverse_delete_non_leaf(bp, node->_left);
 traverse_delete_non_leaf(bp, node->_right);
 zp_pool_release(&bp->_node_allocator, n);
}








    

int zp_broadphase2d_init(zp_broadphase2d *zp_restrict const bp, const float growth_base, float aabb_margin) {
 if(zp_unlikely(zp_pool_init(&bp->_node_allocator, sizeof(zp_broadphase2d_node), 16, growth_base)))
  return -1;
 if(zp_unlikely(zp_bump_init(&bp->_stack, 32)))
  return -1;
 bp->_root = ZP_POOL_NULL_ID;
 bp->_aabb_margin = aabb_margin;
 return 0;
}


void zp_broadphase2d_destroy(zp_broadphase2d *zp_restrict const bp) {
 zp_pool_destroy(&bp->_node_allocator);
 zp_bump_destroy(&bp->_stack);
}


void zp_broadphase2d_update_element(zp_broadphase2d *zp_restrict const bp, const zp_pool_id id, const zp_aabb2d fit) {
 zp_broadphase2d_node *leaf = (zp_broadphase2d_node*)zp_pool_get(&bp->_node_allocator, id);
 if(contains(leaf->_aabb, fit))
  return;
 remove_leaf(bp, id);
 
 leaf->_aabb = expanded(fit, bp->_aabb_margin);
 leaf->_parent = ZP_POOL_NULL_ID;
 leaf->_left = ZP_POOL_NULL_ID;
 leaf->_right = ZP_POOL_NULL_ID;
 insert_leaf(bp, id);
}


zp_pool_id zp_broadphase2d_insert_element(zp_broadphase2d *zp_restrict const bp, const zp_aabb2d fit, const zp_container_id body_id) {
 zp_pool_id new_leaf = create_leaf(bp, expanded(fit, 2.0f), body_id);
 insert_leaf(bp, new_leaf);
 bp->_leaf_size++;
 return new_leaf;
}


void zp_broadphase2d_remove_element(zp_broadphase2d *zp_restrict const bp, const zp_pool_id id) {
 remove_leaf(bp, id);
 zp_pool_release(&bp->_node_allocator, id);
 bp->_leaf_size--;
}


void zp_broadphase2d_traverse_pairs(zp_broadphase2d *zp_restrict const bp, void (*func)(zp_container_id, zp_container_id, void*), void *ptr) {
 zp_broadphase2d_treetranversalinfo tree_info;
 memset(&tree_info, 0, sizeof(zp_broadphase2d_treetranversalinfo));
 traverse_self_recursive(bp, bp->_root, func, ptr, 0, &tree_info);

 float tolerance_factor = 2.5f;
 float log_n = zp_log2((float)bp->_leaf_size);
 uint32_t max_depth = zp_round(tolerance_factor * log_n);
 
 int ca = tree_info._max_depth > max_depth;
 int cb = tree_info._node_visits > (bp->_leaf_size * 20u);
 int cc = tree_info._pair_test > (bp->_leaf_size * 40u);

 if(ca && cb && cc)
  bp->_frame_tolerance++;
 else
  bp->_frame_tolerance = 0;
}



void zp_broadphase2d_optimize_tree(zp_broadphase2d *zp_restrict const bp, void *zp_restrict w) {
 /* recursive delete */
 traverse_delete_non_leaf(bp, bp->_root);
 bp->_root = ZP_POOL_NULL_ID;

 zp_world2d *const world = (zp_world2d*)w;
 zp_body2d *bodies = (zp_body2d*)world->_body_container._bytes;
 const size_t body_size = (size_t)world->_body_container._size;

 if(body_size == 0) {
  return;
 }
  /* query all aabb from the body */
 size_t offset = zp_bump_acquire(&bp->_stack, sizeof(zp_pool_id) * body_size, ZP_MEMORY_ALIGNMENT);
 zp_pool_id *leafs = zp_bump_get(&bp->_stack, offset);
 for(size_t i = 0; i < body_size; i++) {
  leafs[i] = bodies[i]._head._aabb_node;
 }
 
 size_t count = body_size;
 while (count > 1) {

  zp_pool_id parent_id = zp_pool_acquire(&bp->_node_allocator);
  zp_broadphase2d_node *parent = (zp_broadphase2d_node*)zp_pool_get(&bp->_node_allocator, parent_id);
  
  float min_cost = zp_inf();
  int min_i = -1;
  int min_j = -1;

  for(size_t i = 0; i < count; ++i) {
   zp_broadphase2d_node *a = (zp_broadphase2d_node*)zp_pool_get(&bp->_node_allocator, leafs[i]);
   for(size_t j = i + 1; j < count; ++j) {
    zp_broadphase2d_node *b = (zp_broadphase2d_node*)zp_pool_get(&bp->_node_allocator, leafs[j]);
    zp_aabb2d combined = combine(a->_aabb, b->_aabb);
    float cost = perimeter(combined);
    if(cost < min_cost) {
     min_cost = cost;
     min_i = i;
     min_j = j;
    }
   }
  }
  zp_pool_id id_a = leafs[min_i];
  zp_pool_id id_b = leafs[min_j];
 
  zp_broadphase2d_node *node_a = (zp_broadphase2d_node*)zp_pool_get(&bp->_node_allocator, id_a);
  zp_broadphase2d_node *node_b = (zp_broadphase2d_node*)zp_pool_get(&bp->_node_allocator, id_b);
  
  parent->_id = parent_id;
  parent->_left = id_a;
  parent->_right = id_b;
  parent->_parent = ZP_POOL_NULL_ID;
  parent->_aabb = combine(node_a->_aabb, node_b->_aabb);
  
  node_a->_parent = parent_id;
  node_b->_parent = parent_id;

  leafs[min_i] = parent_id;
  leafs[min_j] = leafs[count - 1];
  count--;
 }
 bp->_root = leafs[0];
 zp_bump_reset(&bp->_stack, 0);
}


int zp_broadphase2d_should_rebuild(zp_broadphase2d *zp_restrict const bp, uint8_t frame_tolerance) {
 return bp->_frame_tolerance >= frame_tolerance;
}


