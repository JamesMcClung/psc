
#pragma once

#include "injector_simple.hxx"
#include "const_accessor_simple.hxx"
#include "particles.hxx"
#include "particle_simple.hxx"
#include "particle_indexer.hxx"

#include <iterator>

template<typename Mparticles>
struct MparticlesPatchSimple;

// ======================================================================
// Mparticles

template<typename P>
struct MparticlesSimple : MparticlesBase
{
  using Particle = P;
  using real_t = typename Particle::real_t;
  using Real3 = Vec3<real_t>;
  using BndpParticle = P;
  using Accessor = AccessorSimple<MparticlesSimple>;

  using buf_t = std::vector<Particle>;

  struct Patch
  {
    using iterator = typename buf_t::iterator;
    using const_iterator = typename buf_t::const_iterator;
    
    Patch(MparticlesSimple* mprts, int p)
      : mprts_(mprts),
	p_(p)
      {}
    
    Patch(const Patch&) = delete;
    Patch(Patch&&) = default;

    buf_t&       buf()       { return mprts_->bufs_[p_]; }
    const buf_t& buf() const { return mprts_->bufs_[p_]; }
    
    Particle& operator[](int n) { return buf()[n]; }
    const Particle& operator[](int n) const { return buf()[n]; }
    const_iterator begin() const { return buf().begin(); }
    iterator begin() { return buf().begin(); }
    const_iterator end() const { return buf().end(); }
    iterator end() { return buf().end(); }
    unsigned int size() const { return buf().size(); }
    
    void push_back(const Particle& new_prt)
    {
      // need to copy because we modify it
      auto prt = new_prt;
      checkInPatchMod(prt);
      validCellIndex(prt);
      buf().push_back(prt);
    }
    
    void check() const
    {
      for (auto& prt : buf()) {
	validCellIndex(prt);
      }
    }
    
    // ParticleIndexer functionality
    int cellPosition(real_t xi, int d) const { return mprts_->pi_.cellPosition(xi, d); }
    int validCellIndex(const Particle& prt) const { return mprts_->pi_.validCellIndex(prt.x()); }
    
    void checkInPatchMod(Particle& prt) const { return mprts_->pi_.checkInPatchMod(prt.x()); }
    
    const Grid_t& grid() const { return mprts_->grid(); }
    const MparticlesSimple& mprts() const { return *mprts_; }
    int p() const { return p_; }
    buf_t& bndBuffer() { return buf(); }
    
  public: // FIXME
    MparticlesSimple* mprts_;
  private:
    int p_;
  };

  MparticlesSimple(const Grid_t& grid)
    : MparticlesBase(grid),
      pi_(grid)
  {
    bufs_.resize(grid.n_patches());
  }

  MparticlesSimple(const MparticlesSimple&) = delete;
  MparticlesSimple(MparticlesSimple&& o) = default;

  MparticlesSimple& operator=(MparticlesSimple&& o) = default;

  void reset(const Grid_t& grid) override
  {
    MparticlesBase::reset(grid);
    bufs_ = std::vector<buf_t>(grid.n_patches());
  }

  Patch operator[](int p) const { return {const_cast<MparticlesSimple*>(this), p}; } // FIXME, isn't actually const

  void reserve_all(const std::vector<uint> &n_prts_by_patch)
  {
    for (int p = 0; p < n_patches(); p++) {
      bufs_[p].reserve(n_prts_by_patch[p]);
    }
  }

  void resize_all(const std::vector<uint>& n_prts_by_patch)
  {
    for (int p = 0; p < n_patches(); p++) {
      assert(n_prts_by_patch[p] <= bufs_[p].capacity());
      bufs_[p].resize(n_prts_by_patch[p]);
    }
  }

  void reset() // FIXME, "reset" is used for two very different functions
  {
    for (int p = 0; p < n_patches(); p++) {
      bufs_[p].resize(0);
    }
  }

  std::vector<uint> get_size_all() const override
  {
    std::vector<uint> n_prts_by_patch(n_patches());
    for (int p = 0; p < n_patches(); p++) {
      n_prts_by_patch[p] = bufs_[p].size();
    }
    return n_prts_by_patch;
  }

  int get_n_prts() const override
  {
    int n_prts = 0;
    for (auto const& buf : bufs_) {
      n_prts += buf.size();
    }
    return n_prts;
  }

  const ParticleIndexer<real_t>& particleIndexer() const { return pi_; }
  
  InjectorSimple<MparticlesSimple> injector() { return {*this}; }
  ConstAccessorSimple<MparticlesSimple> accessor() { return {*this}; }
  Accessor accessor_() { return {*this}; }
  
  void check() const
  {
    for (int p = 0; p < n_patches(); p++) {
      (*this)[p].check();
    }
  }

  void dump(const std::string& filename)
  {
    FILE* file = fopen(filename.c_str(), "w");
    assert(file);

    for (int p = 0; p < n_patches(); p++) {
      auto& prts = (*this)[p];
      fprintf(file, "mparticles_dump: p%d n_prts = %d\n", p, prts.size());
      for (int n = 0; n < prts.size(); n++) {
	auto& prt = prts[n];
	fprintf(file, "mparticles_dump: [%d] %g %g %g // %d // %g %g %g // %g\n",
		n, prt.xi, prt.yi, prt.zi, prt.kind_,
		prt.pxi, prt.pyi, prt.pzi, prt.qni_wni_);
      }
    }
    fclose(file);
  }
  
  void define_species(const char *name, double q, double m,
		      double max_local_np, double max_local_nm,
		      double sort_interval, double sort_out_of_place)
  {}
  
  static const Convert convert_to_, convert_from_;
  const Convert& convert_to() override { return convert_to_; }
  const Convert& convert_from() override { return convert_from_; }

private:
  std::vector<buf_t> bufs_;
public: // FIXME
  ParticleIndexer<real_t> pi_;
};

