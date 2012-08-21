#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;
modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," CaDynamics_E2.mod");
    fprintf(stderr," Ca_HVA.mod");
    fprintf(stderr," Ca_LVAst.mod");
    fprintf(stderr," Ih.mod");
    fprintf(stderr," Im.mod");
    fprintf(stderr," K_Pst.mod");
    fprintf(stderr," K_Tst.mod");
    fprintf(stderr," NaTa_t.mod");
    fprintf(stderr," NaTs2_t.mod");
    fprintf(stderr," Nap_Et2.mod");
    fprintf(stderr," SK_E2.mod");
    fprintf(stderr," SKv3_1.mod");
    fprintf(stderr," ampa.mod");
    fprintf(stderr," epsp.mod");
    fprintf(stderr," glutamate.mod");
    fprintf(stderr, "\n");
  }
  _CaDynamics_E2_reg();
  _Ca_HVA_reg();
  _Ca_LVAst_reg();
  _Ih_reg();
  _Im_reg();
  _K_Pst_reg();
  _K_Tst_reg();
  _NaTa_t_reg();
  _NaTs2_t_reg();
  _Nap_Et2_reg();
  _SK_E2_reg();
  _SKv3_1_reg();
  _ampa_reg();
  _epsp_reg();
  _glutamate_reg();
}
