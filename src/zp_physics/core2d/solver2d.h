#ifndef ZP_SOLVER2D_H
#define ZP_SOLVER2D_H

#include "zp_physics/hint.h"
#include "zp_physics/core2d/manifold2d.h"

typedef struct {
 float   _dt;
 float   _inv_dt;
 /* soft constraint values */
 float   _bias_ratio;
 float   _mass_coeff;
 float   _impulse_coeff;
 
} zp_solver_input2d;

void zp_manifold2d_soft_prepare_contact(zp_manifold2d *const zp_restrict m, void *const zp_restrict w, const zp_solver_input2d *const input);
void zp_manifold2d_soft_presolve_contact(zp_manifold2d *const zp_restrict m, void *const zp_restrict w, const zp_solver_input2d *const input);
void zp_manifold2d_soft_solve_contact(zp_manifold2d *const zp_restrict m, void *const zp_restrict w, const zp_solver_input2d *const input);
void zp_manifold2d_soft_relaxation(zp_manifold2d *const zp_restrict m, void *const zp_restrict w, const zp_solver_input2d *const input);


#endif

