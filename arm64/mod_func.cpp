#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;
#if defined(__cplusplus)
extern "C" {
#endif

extern void _model_ion_channels_reg(void);

void modl_reg() {
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");
    fprintf(stderr, " \"model_ion_channels.mod\"");
    fprintf(stderr, "\n");
  }
  _model_ion_channels_reg();
}

#if defined(__cplusplus)
}
#endif
