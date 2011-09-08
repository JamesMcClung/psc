
#include <psc_cuda.h>

EXTERN_C void sort_pairs_host(unsigned int *keys, unsigned int *vals, int n);

// FIXME, specific to 1x8x8, should be in ! .cu, so that cell_map works

static void
sort_find_cell_indices_host(particles_cuda_t *pp, struct psc_patch *patch,
			    unsigned int *d_cnis, unsigned int *d_ids)
{
  int n_part = pp->n_part;
  particles_cuda_dev_t *h_part = &pp->h_part;
  particles_cuda_dev_t *d_part = &pp->d_part;
  unsigned int *h_cnis = (unsigned int *) malloc(n_part * sizeof(*h_cnis));
  unsigned int *h_ids = (unsigned int *) malloc(n_part * sizeof(*h_ids));

  check(cudaMemcpy(h_part->xi4, d_part->xi4, n_part * sizeof(float4),
		   cudaMemcpyDeviceToHost));

  for (int i = 0; i < n_part; i++) {
    particle_cuda_real_t dxi = 1.f / ppsc->dx[0];
    particle_cuda_real_t dyi = 1.f / ppsc->dx[1];
    particle_cuda_real_t dzi = 1.f / ppsc->dx[2];
    particle_cuda_real_t xi[3] = {
      (h_part->xi4[i].x - patch->xb[0]) * dxi,
      (h_part->xi4[i].y - patch->xb[1]) * dyi,
      (h_part->xi4[i].z - patch->xb[2]) * dzi };
    int pos[3];
    for (int d = 0; d < 3; d++) {
      pos[d] = particle_base_real_nint(xi[d]);
    }
    
    int idx = (((pos[2] / 8) * (patch->ldims[1] / 8) + (pos[1] / 8)) << 6) |
      ((pos[2] & 4) << 3) |
      ((pos[2] & 2) << 2) |
      ((pos[2] & 1) << 1) |
      ((pos[1] & 4) << 2) |
      ((pos[1] & 2) << 1) |
      ((pos[1] & 1) << 0);
    h_cnis[i] = idx;
    h_ids[i] = i;
  }

  check(cudaMemcpy(d_cnis, h_cnis, n_part * sizeof(*h_cnis),
		   cudaMemcpyHostToDevice));
  check(cudaMemcpy(d_ids, h_ids, n_part * sizeof(*h_ids),
		   cudaMemcpyHostToDevice));

  free(h_cnis);
  free(h_ids);
}

static void
sort_reorder_host(particles_cuda_t *pp, unsigned int *d_ids)
{
  int n_part = pp->n_part;
  particles_cuda_dev_t *h_part = &pp->h_part;
  particles_cuda_dev_t *d_part = &pp->d_part;
  unsigned int *h_ids = (unsigned int *) malloc(n_part);

  check(cudaMemcpy(h_part->xi4, d_part->xi4, n_part * sizeof(float4),
		   cudaMemcpyDeviceToHost));
  check(cudaMemcpy(h_part->pxi4, d_part->pxi4, n_part * sizeof(float4),
		   cudaMemcpyDeviceToHost));
  check(cudaMemcpy(h_ids, d_ids, n_part * sizeof(*h_ids),
		   cudaMemcpyDeviceToHost));

  // move into new position
  float4 *xi4 = (float4 *) malloc(n_part * sizeof(*xi4));
  float4 *pxi4 = (float4 *) malloc(n_part * sizeof(*pxi4));
  for (int i = 0; i < n_part; i++) {
    xi4[i] = pp->h_part.xi4[h_ids[i]];
    pxi4[i] = pp->h_part.pxi4[h_ids[i]];
  }
  // back to in-place
  memcpy(pp->h_part.xi4, xi4, n_part * sizeof(*xi4));
  memcpy(pp->h_part.pxi4, pxi4, n_part * sizeof(*pxi4));
  
  free(xi4);
  free(pxi4);

  check(cudaMemcpy(d_part->xi4, h_part->xi4, n_part * sizeof(float4),
		   cudaMemcpyHostToDevice));
  check(cudaMemcpy(d_part->pxi4, h_part->pxi4, n_part * sizeof(float4),
		   cudaMemcpyHostToDevice));
  free(h_ids);
}

EXTERN_C void
sort_patch(int p, particles_cuda_t *pp)
{
  struct psc_patch *patch = &ppsc->patch[p];

  unsigned int *d_cnis, *d_ids;
  check(cudaMalloc((void **) &d_cnis, pp->n_part * sizeof(*d_cnis)));
  check(cudaMalloc((void **) &d_ids, pp->n_part * sizeof(*d_ids)));

  sort_find_cell_indices_host(pp, patch, d_cnis, d_ids);
  sort_pairs_host(d_cnis, d_ids, pp->n_part);
  sort_reorder_host(pp, d_ids);

  check(cudaFree(d_cnis));
  check(cudaFree(d_ids));
}
