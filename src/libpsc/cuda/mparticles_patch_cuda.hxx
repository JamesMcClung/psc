
#pragma once

// FIXME?  we have two pretty similar versions of ConstAccessorCuda here,
// and probably only one should survive.  one version copies from
// device "on demand", whereas the other one copies a whole patch
// worth of data just in case.  Obviously, there can be use cases for
// either, but this needs some thinking as to what's actually needed
// for this code.  The on-demand version might be useful to serve as a
// template for a modifiable accesser, if we ever want to go there.

// ======================================================================
// ConstAccessorCuda

template<typename Mparticles>
struct ConstAccessorCuda
{
  using Particle = typename Mparticles::Particle;
  using real_t = typename Particle::real_t;
  using Real3 = typename Particle::Real3;

  struct const_accessor
  {
    using Double3 = Vec3<double>;
    
    const_accessor(const Particle& prt, const Mparticles& mprts, int p)
      : prt_{prt}, mprts_{mprts}, p_{p}
    {}

    Real3 x()   const { return prt_.x(); }
    Real3 u()   const { return prt_.u(); }
    real_t w()  const { return prt_.qni_wni() / mprts_.grid().kinds[prt_.kind()].q; }
    real_t qni_wni() const { return prt_.qni_wni(); }
    int kind()  const { return prt_.kind(); }
    
    Double3 position() const
    {
      auto& patch = mprts_.grid().patches[p_];
      
      return patch.xb + Double3(prt_.x());
    }

    operator const Particle& () const { return prt_; }
    
  private:
    const Particle prt_;
    const Mparticles& mprts_;
    const int p_;
  };
  
  struct Patch
  {
    struct const_iterator : std::iterator<std::random_access_iterator_tag,
      const_accessor,  // value type
      ptrdiff_t,       // difference type
      const_accessor*, // pointer type
      const_accessor&> // reference type
    
    {
      const_iterator(const Patch& patch, uint n)
	: patch_{patch}, n_{n}
      {}

      bool operator==(const_iterator other) const { return n_ == other.n_; }
      bool operator!=(const_iterator other) const { return !(*this == other); }
      
      const_iterator& operator++() { n_++; return *this; }
      const_iterator operator++(int) { auto retval = *this; ++(*this); return retval; }
      const_accessor operator*() { return {patch_.data_[n_], patch_.mprts_, patch_.p_}; }
      
    private:
      const Patch patch_;
      uint n_;
    };
    
    Patch(const Mparticles& mprts, int p)
      : mprts_{mprts}, p_{p}, data_{const_cast<Mparticles&>(mprts).get_particles(p)}
    // FIXME, const hacking around reorder may change state...
    {}

    // FIXME, implicit copy ctor copies entire array, and that happens all the time when
    // making a const_iterator, which is rather bad
      
    const_iterator begin() const { return {*this, 0}; }
    const_iterator end()   const { return {*this, uint(data_.size())}; }
    const_accessor operator[](int n) const { return {data_[n], mprts_, p_}; }
    uint size() const { return data_.size(); }
    
  private:
    const Mparticles& mprts_;
    int p_;
    const std::vector<Particle> data_;
  };

  ConstAccessorCuda(Mparticles& mprts)
    : mprts_{mprts}
  {}

  Patch operator[](int p) { return {mprts_, p}; }

private:
  Mparticles& mprts_;
};

// ======================================================================
// ConstParticleAccessorCuda_

template<typename Mparticles>
struct ConstParticleAccessorCuda_
{
  using Particle = typename Mparticles::Particle;
  using real_t = typename Particle::real_t;
  using Real3 = Vec3<real_t>;
  using Double3 = Vec3<double>;

  ConstParticleAccessorCuda_(const Particle& prt, const Mparticles& mprts, int p)
    : prt_{prt}, mprts_{mprts}, p_{p}
  {}
  
  Real3 x()   const { return prt_.x(); }
  Real3 u()   const { return prt_.u(); }
  real_t w()  const { return prt_.qni_wni() / mprts_.grid().kinds[prt_.kind()].q; }
  real_t qni_wni() const { return prt_.qni_wni(); }
  int kind()  const { return prt_.kind(); }
  
  Double3 position() const
  {
    auto& patch = mprts_.grid().patches[p_];
    
    return patch.xb + Double3(prt_.x());
  }
  
  operator const Particle& () const { return prt_; }
  
private:
  Particle prt_;
  const Mparticles& mprts_;
  const int p_;
};
  

// ======================================================================
// ConstAccessorCuda_

template<typename Mparticles>
struct ConstAccessorCuda_
{
  using Particle = typename Mparticles::Particle;
  using real_t = typename Particle::real_t;
  using Real3 = typename Particle::Real3;

  using const_accessor = ConstParticleAccessorCuda_<Mparticles>;
  
  struct Patch
  {
    struct const_iterator : std::iterator<std::random_access_iterator_tag,
					  const_accessor,  // value type
					  ptrdiff_t,       // difference type
					  const_accessor*, // pointer type
					  const_accessor&> // reference type
      
    {
      const_iterator(const Mparticles& mprts, int p, uint n)
	: mprts_{mprts}, p_{p}, n_{n}
      {}
	
      bool operator==(const_iterator other) const { return n_ == other.n_; }
      bool operator!=(const_iterator other) const { return !(*this == other); }
	
      const_iterator& operator++() { n_++; return *this; }
      const_iterator operator++(int) { auto retval = *this; ++(*this); return retval; }
      const_accessor operator*() { return {const_cast<Mparticles&>(mprts_).get_particle(p_, n_), mprts_, p_}; } // FIXME constness
	
    private:
      const Mparticles& mprts_;
      const int p_;
      uint n_;
    };
    
    Patch(const Mparticles& mprts, int p)
      : mprts_{mprts}, p_{p}
    {}

    const_iterator begin() const { return {mprts_, p_, 0}; }
    const_iterator end()   const { return {mprts_, p_, size()}; };
    const_accessor operator[](int n) const { return {const_cast<Mparticles&>(mprts_).get_particle(p_, n), mprts_, p_}; }

    uint size() const
    {
      auto n_prts_by_patch = const_cast<Mparticles&>(mprts_).get_size_all();
      return n_prts_by_patch[p_];
    }

      
  private:
    const Mparticles& mprts_;
    const int p_;
  };

  ConstAccessorCuda_(Mparticles& mprts)
    : mprts_{mprts}
  {}

  Patch operator[](int p) { return {mprts_, p}; }

private:
  Mparticles& mprts_;
};

// ======================================================================
// ConstPatchCuda

template<typename _Mparticles>
struct ConstPatchCuda
{
  using Mparticles = _Mparticles;

  ConstPatchCuda(const Mparticles& mprts, int p)
    : mprts_{mprts}, p_(p)
  {}
  
private:
  const Mparticles& mprts_;
  int p_;
};

