#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;
#if defined(__cplusplus)
extern "C" {
#endif

extern void _km_reg(void);
extern void _model_ion_channels_reg(void);
extern void _v_gaba_b_reg(void);
extern void _v_nmda_reg(void);

void modl_reg() {
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");
    fprintf(stderr, " \"km.mod\"");
    fprintf(stderr, " \"model_ion_channels.mod\"");
    fprintf(stderr, " \"v_gaba_b.mod\"");
    fprintf(stderr, " \"v_nmda.mod\"");
    fprintf(stderr, "\n");
  }
  _km_reg();
  _model_ion_channels_reg();
  _v_gaba_b_reg();
  _v_nmda_reg();
}

#if defined(__cplusplus)
}
#endif
