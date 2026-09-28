#include "zp_physics/core2d/island2d.h"
#include <assert.h>
#include <string.h>

/**********************************************
*               ISLAND OPERATION              *
***********************************************/


zp_noinline size_t find(zp_island2d *zp_restrict const island, size_t x_id) {
 size_t mx = x_id;
 zp_island2d_id *x = (zp_island2d_id*)zp_bump_get(&island->_stack, mx);
 while(x->_parent != mx) {
  zp_island2d_id *x_parent = (zp_island2d_id*)zp_bump_get(&island->_stack, x->_parent);
  x->_parent = x_parent->_parent;
  mx = x_parent->_allocation;
  x = (zp_island2d_id*)zp_bump_get(&island->_stack, mx);
 }
 return mx;
}




void zp_island2d_init(zp_island2d *zp_restrict const island, size_t reserve, float growth_base) {
 zp_container_init(&island->_island_ids, sizeof(zp_island_representative), reserve, growth_base);
 zp_bump_init(&island->_stack, 512, growth_base);
}


void zp_island2d_destroy(zp_island2d *zp_restrict const island) {
 zp_container_destroy(&island->_island_ids);
 zp_bump_destroy(&island->_stack);
}


void zp_island2d_merge(zp_island2d *zp_restrict const island, zp_container *zp_restrict const container, zp_body2d *zp_restrict const body_a, zp_body2d *zp_restrict const body_b) {
 (void)island;
 (void)container;
 (void)body_a;
 (void)body_b;
 /*
 IslandId* ia = find(a->island);
 IslandId* ib = find(b->island);

    if (ia == ib)
        return ia;

    // Larger island becomes parent.
    if (ia->size < ib->size)
        std::swap(ia, ib);

    // ------------------------------------------------
    // DSU merge
    // ------------------------------------------------

    ib->parent = ia;
    ia->size += ib->size;
    
    // ------------------------------------------------
    // O(1) linked-list splice
    // ------------------------------------------------

    ia->tail->next = ib->head;
    ib->head->prev = ia->tail;

    ia->tail = ib->tail;

    return ia;*/
}

    
void zp_island2d_reset(zp_island2d *zp_restrict const island) {
 zp_bump_reset(&island->_stack, 0);
}

/*
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>


typedef struct Body Body;
typedef struct Island Island;
typedef struct Contact Contact;



struct Body
{

    Body *next;


    Island *island;


    float x;
    float y;

    float vx;
    float vy;

    float mass;
};


struct Island
{

    Island *parent;

    Island *prev_root;
    Island *next_root;

    Body *head;
    Body *tail;

    size_t body_count;


    int sleeping;
};


typedef struct
{
    Island *root_list;

    size_t island_count;
} IslandWorld;


static void island_world_init(IslandWorld *world)
{
    world->root_list = NULL;
    world->island_count = 0;
}



static void island_root_add(IslandWorld *world, Island *island)
{
    island->parent = NULL;

    island->prev_root = NULL;
    island->next_root = world->root_list;

    if (world->root_list != NULL)
        world->root_list->prev_root = island;

    world->root_list = island;

    world->island_count++;
}



static void island_root_remove(IslandWorld *world, Island *island)
{
    if (island->prev_root != NULL)
        island->prev_root->next_root = island->next_root;
    else
        world->root_list = island->next_root;

    if (island->next_root != NULL)
        island->next_root->prev_root = island->prev_root;

    island->prev_root = NULL;
    island->next_root = NULL;

    world->island_count--;
}



static Island *island_create(IslandWorld *world, Body *body)
{
    Island *island = (Island*)malloc(sizeof(*island));

    if (island == NULL)
        return NULL;

    island->parent = NULL;

    island->prev_root = NULL;
    island->next_root = NULL;

    island->head = body;
    island->tail = body;

    island->body_count = 1;
    island->sleeping = 0;

    body->next = NULL;
    body->island = island;

    island_root_add(world, island);

    return island;
}



static Island *island_find(Island *island)
{
    Island *root;

    if (island == NULL)
        return NULL;


    root = island;

    while (root->parent != NULL)
        root = root->parent;


    while (island->parent != NULL)
    {
        Island *next = island->parent;

        island->parent = root;
        island = next;
    }

    return root;
}


static Island *island_merge(
    IslandWorld *world,
    Island *a,
    Island *b)
{
    if (a == NULL || b == NULL)
        return a ? a : b;


    if (a == b)
        return a;


    if (a->parent != NULL || b->parent != NULL)
        return NULL;

    if (b->tail != NULL)
    {
        b->tail->next = a->head;
        b->tail = a->tail;
    }
    else
    {
        b->head = a->head;
        b->tail = a->tail;
    }

    b->body_count += a->body_count;


    a->parent = b;

    island_root_remove(world, a);

    return b;
}



static Island *island_connect(
    IslandWorld *world,
    Body *a,
    Body *b)
{
    Island *ia;
    Island *ib;

    if (a == NULL || b == NULL)
        return NULL;

    ia = island_find(a->island);
    ib = island_find(b->island);


    if (ia == ib)
        return ia;

    return island_merge(world, ia, ib);
}



struct Contact
{
    Body *body_a;
    Body *body_b;

    Contact *next;
};



static void island_process_contacts(
    IslandWorld *world,
    Contact *contacts)
{
    Contact *contact;

    for (contact = contacts;
         contact != NULL;
         contact = contact->next)
    {
        island_connect(
            world,
            contact->body_a,
            contact->body_b);
    }
}



static void integrate_island(Island *island, float dt)
{
    Body *body;



    for (body = island->head;
         body != NULL;
         body = body->next)
    {
        body->x += body->vx * dt;
        body->y += body->vy * dt;
    }
}



static void integrate_all_islands(
    IslandWorld *world,
    float dt)
{
    Island *island;

    for (island = world->root_list;
         island != NULL;
         island = island->next_root)
    {
        integrate_island(island, dt);
    }
}



static void print_islands(const IslandWorld *world)
{
    const Island *island;
    size_t island_index = 0;

    printf("Islands: %zu\n", world->island_count);

    for (island = world->root_list;
         island != NULL;
         island = island->next_root)
    {
        const Body *body;

        printf("  Island %zu (%zu bodies): ",
               island_index++,
               island->body_count);

        for (body = island->head;
             body != NULL;
             body = body->next)
        {
            printf("%p ", (void *)body);
        }

        printf("\n");
    }
}



int main(void)
{
    IslandWorld world;

    Body a = {0};
    Body b = {0};
    Body c = {0};
    Body d = {0};

    Contact c1;
    Contact c2;
    Contact c3;

    island_world_init(&world);

 

    if (!island_create(&world, &a) ||
        !island_create(&world, &b) ||
        !island_create(&world, &c) ||
        !island_create(&world, &d))
    {
        return 1;
    }



    c1.body_a = &a;
    c1.body_b = &b;
    c1.next = &c2;

    c2.body_a = &b;
    c2.body_b = &c;
    c2.next = &c3;

    c3.body_a = &a;
    c3.body_b = &c;
    c3.next = NULL;

    island_process_contacts(&world, &c1);

    print_islands(&world);


    integrate_all_islands(&world, 1.0f);

    return 0;
}


*/