
#ifndef PSC_PARTICLES_AS_CBE_H
#define PSC_PARTICLES_AS_CBE_H

#include "psc_particles_cbe.h"

typedef particle_cbe_t particle_t;
typedef particles_cbe_t particles_t;
typedef mparticles_cbe_t mparticles_t;

#define particles_get          particles_cbe_get
#define particles_put          particles_cbe_put
#define particles_get_one      particles_cbe_get_one

#endif


/// \file psc_particles_as_cbe.h Some defines for cell be particles.
///
/// This file is strictly temporary. If I can get a couple issues resolved, 
/// (one will involving forcing alignment in the particle arrays) then 
/// I can just use the c particles instead of having a unique particle type.