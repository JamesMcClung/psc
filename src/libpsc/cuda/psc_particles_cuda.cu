
#include "mparticles_cuda.hxx"
#include "psc.h"
#include "cuda_mparticles.cuh"
#include "psc_particles_single.h"
#include "psc_particles_double.h"

#include <mrc_io.h>

#if 0
#define dprintf(...) mprintf(__VA_ARGS__)
#else
#define dprintf(...) do {} while (0)
#endif

// ======================================================================
// MparticleCuda implementation

template<typename BS>
MparticlesCuda<BS>::MparticlesCuda(const Grid_t& grid)
  : MparticlesBase(grid),
    pi_(grid)
{
  dprintf("CMPRTS: ctor\n");
  cmprts_ = new cuda_mparticles<BS>(grid);
}

template<typename BS>
MparticlesCuda<BS>::~MparticlesCuda()
{
  dprintf("CMPRTS: dtor\n");
  delete cmprts_;
}

template<typename BS>
std::vector<uint> MparticlesCuda<BS>::get_size_all() const
{
  dprintf("CMPRTS: get_size_all\n");
  auto n_prts_by_patch = cmprts_->get_size_all();
  for (int p = 0; p < cmprts_->n_patches(); p++) {
    dprintf("  p %d: %d\n", p, n_prts_by_patch[p]);
  }
  return n_prts_by_patch;
}

template<typename BS>
int MparticlesCuda<BS>::get_n_prts() const
{
  dprintf("CMPRTS: get_n_prts\n");
  return cmprts_->get_n_prts();
}

template<typename BS>
void MparticlesCuda<BS>::reset(const Grid_t& grid)
{
  this->~MparticlesCuda();
  new(this) MparticlesCuda(grid);
}

template<typename BS>
void MparticlesCuda<BS>::inject(const std::vector<particle_t>& buf, const std::vector<uint>& buf_n_by_patch)
{
  dprintf("CMPRTS: inject\n");
  cmprts_->inject(buf, buf_n_by_patch);
}

template<typename BS>
std::vector<typename MparticlesCuda<BS>::particle_t> MparticlesCuda<BS>::get_particles(int p) const
{
  dprintf("CMPRTS: get_particles\n");
  return cmprts_->get_particles(p);
}

template<typename BS>
typename MparticlesCuda<BS>::particle_t MparticlesCuda<BS>::get_particle(int p, int n) const
{
  dprintf("CMPRTS: get_particle\n");
  return cmprts_->get_particle(p, n);
}

template<typename BS>
void MparticlesCuda<BS>::dump(const std::string& filename)
{
  cmprts_->dump(filename);
}

template<typename BS>
bool MparticlesCuda<BS>::check_after_push()
{
  return cmprts_->check_bidx_after_push();
}

// ======================================================================
// conversion

template<typename MparticlesCuda, typename MP>
static void copy_from(MparticlesBase& mprts_base, MparticlesBase& mprts_other_base)
{
  auto& mprts = dynamic_cast<MparticlesCuda&>(mprts_base);
  auto& mprts_other = dynamic_cast<MP&>(mprts_other_base);
  auto n_prts_by_patch = mprts_other.get_size_all();
  //mp.reserve_all(n_prts_by_patch); FIXME, would still be a good hint for the injector

  auto accessor = mprts_other.accessor();
  auto inj = mprts.injector();
  for (int p = 0; p < mprts.n_patches(); p++) {
    auto injector = inj[p];
    for (auto prt: accessor[p]) {
      using real_t = typename MparticlesCuda::real_t;
      using Real3 = typename MparticlesCuda::Real3;
      injector.raw({Real3(prt.x()), Real3(prt.u()), real_t(prt.qni_wni()), prt.kind()});
    }
  }
}

template<typename MparticlesCuda, typename MP>
static void copy_to(MparticlesBase& mprts_base, MparticlesBase& mprts_other_base)
{
  auto& mprts = dynamic_cast<MparticlesCuda&>(mprts_base);
  auto& mprts_other = dynamic_cast<MP&>(mprts_other_base);
  auto n_prts_by_patch = mprts.get_size_all();
  mprts_other.reserve_all(n_prts_by_patch);
  mprts_other.resize_all(n_prts_by_patch);

  auto accessor = mprts.accessor(); // FIXME, should we use this in the first place?
  for (int p = 0; p < mprts.n_patches(); p++) {
    int n = 0;
    for (auto prt: accessor[p]) {
      using real_t = typename MP::real_t;
      using Real3 = typename MP::Real3;
      mprts_other[p][n++] = {Real3(prt.x()), Real3(prt.u()), real_t(prt.qni_wni()), prt.kind()};
    }
  }
}

// ======================================================================
// conversion to "single"/"double"

template<typename BS>
const MparticlesCuda<BS>::Convert MparticlesCuda<BS>::convert_to_ = {
  { std::type_index(typeid(MparticlesSingle)), copy_to<MparticlesCuda<BS>, MparticlesSingle>   },
  { std::type_index(typeid(MparticlesDouble)), copy_to<MparticlesCuda<BS>, MparticlesDouble>   },
};

template<typename BS>
const MparticlesCuda<BS>::Convert MparticlesCuda<BS>::convert_from_ = {
  { std::type_index(typeid(MparticlesSingle)), copy_from<MparticlesCuda<BS>, MparticlesSingle>   },
  { std::type_index(typeid(MparticlesDouble)), copy_from<MparticlesCuda<BS>, MparticlesDouble>   },
};

// ======================================================================
// psc_mparticles "cuda"

#if 0
#ifdef HAVE_LIBHDF5_HL

// FIXME. This is a rather bad break of proper layering, HDF5 should be all
// mrc_io business. OTOH, it could be called flexibility...

#include <hdf5.h>
#include <hdf5_hl.h>

#define H5_CHK(ierr) assert(ierr >= 0)
#define CE assert(ierr == 0)

// ----------------------------------------------------------------------
// psc_mparticles_cuda_write

static void
psc_mparticles_cuda_write(struct psc_mparticles *_mprts, struct mrc_io *io)
{
  PscMparticlesCuda mprts(_mprts);
  int ierr;
  
  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);

  auto n_prts_by_patch = mprts->get_size_all();

  hid_t group = H5Gopen(h5_file, mrc_io_obj_path(io, _mprts), H5P_DEFAULT); H5_CHK(group);
  uint off = 0;
  // FIXME, reorder first if necessary
  for (int p = 0; p < mprts.n_patches(); p++) {
    char pname[10];
    sprintf(pname, "p%d", p);
    hid_t pgroup = H5Gcreate(group, pname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); H5_CHK(pgroup);
    int n_prts = n_prts_by_patch[p];
    ierr = H5LTset_attribute_int(pgroup, ".", "n_prts", &n_prts, 1); CE;
    if (n_prts > 0) {
      float4 *xi4  = (float4 *) calloc(n_prts, sizeof(*xi4));
      float4 *pxi4 = (float4 *) calloc(n_prts, sizeof(*pxi4));
      
      mprts->from_device(xi4, pxi4, n_prts, off);
      
      hsize_t hdims[2];
      hdims[0] = n_prts; hdims[1] = 4;
      ierr = H5LTmake_dataset_float(pgroup, "xi4", 2, hdims, (float *) xi4); CE;
      ierr = H5LTmake_dataset_float(pgroup, "pxi4", 2, hdims, (float *) pxi4); CE;
      
      free(xi4);
      free(pxi4);
    }
    ierr = H5Gclose(pgroup); CE;
    off += n_prts;
  }
  ierr = H5Gclose(group); CE;
}

// ----------------------------------------------------------------------
// psc_mparticles_cuda_read

static void
psc_mparticles_cuda_read(struct psc_mparticles *_mprts, struct mrc_io *io)
{
  PscMparticlesCuda mprts(_mprts);

  psc_mparticles_read_super(_mprts, io);

  int ierr;
  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);
  hid_t group = H5Gopen(h5_file, mrc_io_obj_path(io, _mprts), H5P_DEFAULT); H5_CHK(group);

  int n_prts_by_patch[mprts.n_patches()];

  for (int p = 0; p < mprts.n_patches(); p++) {
    char pname[10];
    sprintf(pname, "p%d", p);
    hid_t pgroup = H5Gopen(group, pname, H5P_DEFAULT); H5_CHK(pgroup);
    int n_prts;
    ierr = H5LTget_attribute_int(pgroup, ".", "n_prts", &n_prts); CE;
    n_prts_by_patch[p] = n_prts;
    ierr = H5Gclose(pgroup); CE;
  }

  psc_mparticles_setup(_mprts);
  mprts->reserve_all(n_prts_by_patch);

  uint off = 0;
  for (int p = 0; p < mprts.n_patches(); p++) {
    char pname[10];
    sprintf(pname, "p%d", p);
    hid_t pgroup = H5Gopen(group, pname, H5P_DEFAULT); H5_CHK(pgroup);
    int n_prts = n_prts_by_patch[p];
    if (n_prts > 0) {
      float4 *xi4  = (float4*) calloc(n_prts, sizeof(float4));
      float4 *pxi4 = (float4*) calloc(n_prts, sizeof(float4));
      
      ierr = H5LTread_dataset_float(pgroup, "xi4", (float *) xi4); CE;
      ierr = H5LTread_dataset_float(pgroup, "pxi4", (float *) pxi4); CE;
      
      PscMparticlesCuda mprts = PscMparticlesCuda(_mprts);
      mprts->to_device(xi4, pxi4, n_prts, off);
      
      free(xi4);
      free(pxi4);
    }

    ierr = H5Gclose(pgroup); CE;
    off += n_prts;
  }

  ierr = H5Gclose(group); CE;
  mprts->setup_internals();
}

#endif
#endif

// ----------------------------------------------------------------------
// psc_mparticles: subclass "cuda"
  
template struct MparticlesCuda<BS144>;
template struct MparticlesCuda<BS444>;
