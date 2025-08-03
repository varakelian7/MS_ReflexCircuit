#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;
#if defined(__cplusplus)
extern "C" {
#endif

extern void _gaba_b_reg(void);
extern void _km_reg(void);
extern void _model_ion_channels_reg(void);
extern void _nmda_reg(void);

void modl_reg() {
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");
    fprintf(stderr, " \"gaba_b.mod\"");
    fprintf(stderr, " \"km.mod\"");
    fprintf(stderr, " \"model_ion_channels.mod\"");
    fprintf(stderr, " \"nmda.mod\"");
    fprintf(stderr, "\n");
  }
  _gaba_b_reg();
  _km_reg();
  _model_ion_channels_reg();
  _nmda_reg();
}

#if defined(__cplusplus)
}
#endif
