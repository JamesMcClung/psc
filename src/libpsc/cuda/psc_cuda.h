
#ifndef PSC_CUDA_H
#define PSC_CUDA_H

#include "psc.h"

#include <assert.h>
#include <math.h>
#include <psc.h>

#include <psc_particles.h>

// ======================================================================

void psc_push_particles_cuda_push_yz_a(struct psc_push_particles *push,
				       mparticles_base_t *particles_base,
				       mfields_base_t *flds_base);

void psc_push_particles_cuda_push_yz_b(struct psc_push_particles *push,
				       mparticles_base_t *particles_base,
				       mfields_base_t *flds_base);

void psc_push_particles_cuda_push_yz(struct psc_push_particles *push,
				     mparticles_base_t *particles_base,
				     mfields_base_t *flds_base);

// ======================================================================

#define check(a) do { int ierr = a; if (ierr != cudaSuccess) fprintf(stderr, "IERR = %d (%d)\n", ierr, cudaSuccess); assert(ierr == cudaSuccess); } while(0)

// ======================================================================

#define DECLARE_CUDA(pfx)                                               \
  EXTERN_C void pfx##_set_constants(particles_cuda_t *pp,		\
				    fields_cuda_t *pf);			\
  EXTERN_C void pfx##_cuda_push_part_p1(particles_cuda_t *pp,           \
                                        fields_cuda_t *pf,              \
                                        real **d_scratch);              \
  EXTERN_C void pfx##_cuda_push_part_p2(particles_cuda_t *pp,           \
                                        fields_cuda_t *pf);             \
  EXTERN_C void pfx##_cuda_push_part_p3(particles_cuda_t *pp,           \
                                        fields_cuda_t *pf,              \
                                        real *d_scratch);               \
  EXTERN_C void pfx##_cuda_push_part_p4(particles_cuda_t *pp,           \
                                        fields_cuda_t *pf,              \
                                        real *d_scratch);               \
  EXTERN_C void pfx##_cuda_push_part_p5(particles_cuda_t *pp,           \
                                        fields_cuda_t *pf,              \
                                        real *d_scratch);               \

DECLARE_CUDA(yz);
DECLARE_CUDA(yz2);
DECLARE_CUDA(yz3);

EXTERN_C void yz_a_set_constants(particles_cuda_t *pp, fields_cuda_t *pf);
EXTERN_C void yz_b_set_constants(particles_cuda_t *pp, fields_cuda_t *pf);
EXTERN_C void __cuda_push_part_yz_a(particles_cuda_t *pp, fields_cuda_t *pf);
EXTERN_C void __cuda_push_part_yz_b(particles_cuda_t *pp, fields_cuda_t *pf);
EXTERN_C void __cuda_push_part_yz_b3(particles_cuda_t *pp, fields_cuda_t *pf);

struct d_particle {
  real xi[3];
  real qni_div_mni;
  real pxi[3];
  real qni_wni;
};

#define THREADS_PER_BLOCK 128

#define LOAD_PARTICLE(pp, d_p, n) do {					\
    (pp).xi[0]       = d_p.xi4[n].x;					\
    (pp).xi[1]       = d_p.xi4[n].y;					\
    (pp).xi[2]       = d_p.xi4[n].z;					\
    (pp).qni_div_mni = d_p.xi4[n].w;					\
    (pp).pxi[0]      = d_p.pxi4[n].x;					\
    (pp).pxi[1]      = d_p.pxi4[n].y;					\
    (pp).pxi[2]      = d_p.pxi4[n].z;					\
    (pp).qni_wni     = d_p.pxi4[n].w;					\
} while (0)

#define STORE_PARTICLE_POS(pp, d_p, n) do {				\
    d_p.xi4[n].x = (pp).xi[0];						\
    d_p.xi4[n].y = (pp).xi[1];						\
    d_p.xi4[n].z = (pp).xi[2];						\
    d_p.xi4[n].w = (pp).qni_div_mni;					\
} while (0)

#define STORE_PARTICLE_MOM(pp, d_p, n) do {				\
    d_p.pxi4[n].x = (pp).pxi[0];					\
    d_p.pxi4[n].y = (pp).pxi[1];					\
    d_p.pxi4[n].z = (pp).pxi[2];					\
    d_p.pxi4[n].w = (pp).qni_wni;					\
} while (0)

#endif
