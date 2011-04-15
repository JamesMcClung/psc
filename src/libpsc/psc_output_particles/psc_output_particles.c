
#include "psc_output_particles_private.h"

// ======================================================================
// forward to subclass

void
psc_output_particles_run(struct psc_output_particles *output_particles,
			 mparticles_base_t *particles)
{
  struct psc_output_particles_ops *ops = psc_output_particles_ops(output_particles);
  assert(ops->run);
  ops->run(output_particles, particles);
}

// ======================================================================
// psc_output_particles_init

static void
psc_output_particles_init()
{
  mrc_class_register_subclass(&mrc_class_psc_output_particles,
			      &psc_output_particles_none_ops);
  mrc_class_register_subclass(&mrc_class_psc_output_particles,
			      &psc_output_particles_fortran_ops);
}

// ======================================================================
// psc_output_particles class

struct mrc_class_psc_output_particles mrc_class_psc_output_particles = {
  .name             = "psc_output_particles",
  .size             = sizeof(struct psc_output_particles),
  .init             = psc_output_particles_init,
};
