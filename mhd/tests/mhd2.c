
//#define BOUNDS_CHECK

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_step.h"
#include "ggcm_mhd_diag.h"
#include "ggcm_mhd_ic.h"

#include <mrc_ts.h>
#include <mrc_ts_monitor.h>
#include <mrc_fld.h>
#include <mrc_domain.h>

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h> 
#include <assert.h>

static void
ggcm_mhd_cweno_create(struct ggcm_mhd *mhd)
{
  ggcm_mhd_default_box(mhd);

  /* set defaults for coord arrays */
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  mrc_crds_set_type(crds, "gaussian_2D");
  mrc_crds_set_param_int(crds, "sw", SW_2);   // 'stencil width' 
  mrc_crds_set_param_float3(crds, "l", (float[3]) {  0.0, 0.0, -1.0 });
  mrc_crds_set_param_float3(crds, "h", (float[3]) {  2.*M_PI, 2.*M_PI,  1.0 });
}

static struct ggcm_mhd_ops ggcm_mhd_cweno_ops = {
  .name             = "cweno",
  .create           = ggcm_mhd_cweno_create,
};

// ======================================================================

extern struct ggcm_mhd_diag_ops ggcm_mhd_diag_c_ops;

extern struct ggcm_mhd_ic_ops ggcm_mhd_ic_mirdip_ops;
extern struct ggcm_mhd_ic_ops ggcm_mhd_ic_mirdip2_ops;
extern struct ggcm_mhd_ic_ops ggcm_mhd_ic_mirdip3_ops;
extern struct ggcm_mhd_ic_ops ggcm_mhd_ic_otzi_ops;
extern struct ggcm_mhd_ic_ops ggcm_mhd_ic_ot_ops;
extern struct ggcm_mhd_ic_ops ggcm_mhd_ic_harris_ops;
extern struct ggcm_mhd_ic_ops ggcm_mhd_ic_fadeev_ops;
extern struct ggcm_mhd_ic_ops ggcm_mhd_ic_hydroblast_ops;
extern struct ggcm_mhd_ic_ops ggcm_mhd_ic_mhdblast_ops;
extern struct ggcm_mhd_ic_ops ggcm_mhd_ic_ici_ops;
extern struct ggcm_mhd_ic_ops ggcm_mhd_ic_harris;
extern struct ggcm_mhd_ic_ops ggcm_mhd_ic_wave_sound_ops;
extern struct ggcm_mhd_ic_ops ggcm_mhd_ic_wave_alfven_ops;

int
main(int argc, char **argv)
{
  mrc_class_register_subclass(&mrc_class_ggcm_mhd, &ggcm_mhd_cweno_ops);  
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_diag, &ggcm_mhd_diag_c_ops);

  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_otzi_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_hydroblast_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_mhdblast_ops);    
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_ici_ops); 
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_wave_sound_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_ic, &ggcm_mhd_ic_wave_alfven_ops);
 
  return ggcm_mhd_main(&argc, &argv);
}

