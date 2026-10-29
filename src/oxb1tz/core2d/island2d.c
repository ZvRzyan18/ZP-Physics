#include "oxb1tz/core2d/island2d.h"

/*
 Disjoint Set Union Algorithm
 a fast graph merge with path compression
*/

oxb_pool_id oxb_island2d_find(oxb_island2d *oxb_restrict const island, const oxb_pool_id h) {
 oxb_island2d_node *node = oxb_pool_get(&island->_node_pool, h);
 if(!oxb_pool_id_isequal(node->_leader, h)) 
  node->_leader = oxb_island2d_find(island, node->_leader);
 return node->_leader;
}


void oxb_island2d_merge(oxb_island2d *oxb_restrict const island, const oxb_pool_id a, const oxb_pool_id b) {
 oxb_pool_id ha = oxb_island2d_find(island, a);
 oxb_pool_id hb = oxb_island2d_find(island, b);

 oxb_island2d_node *root_a = oxb_pool_get(&island->_node_pool, ha);
 oxb_island2d_node *root_b = oxb_pool_get(&island->_node_pool, hb);

 if(oxb_pool_id_isequal(ha, hb))
  return;

 if(root_a->_rank < root_b->_rank) {
  root_a->_leader = root_b->_allocation;
 } else if(root_a->_rank > root_b->_rank) {
  root_b->_leader = root_a->_allocation;
 } else {
  root_b->_leader = root_a->_allocation;
  root_a->_rank++;
 }
}



