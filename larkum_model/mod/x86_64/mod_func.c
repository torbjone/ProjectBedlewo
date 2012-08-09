#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;
modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," Kdr.mod");
    fprintf(stderr," SlowCa.mod");
    fprintf(stderr," ampa.mod");
    fprintf(stderr," cad2.mod");
    fprintf(stderr," glutamate.mod");
    fprintf(stderr," h.mod");
    fprintf(stderr," h2.mod");
    fprintf(stderr," hh3.mod");
    fprintf(stderr," ih.mod");
    fprintf(stderr," it2.mod");
    fprintf(stderr," kap.mod");
    fprintf(stderr," kca.mod");
    fprintf(stderr," kdf.mod");
    fprintf(stderr," kdr2.mod");
    fprintf(stderr," km.mod");
    fprintf(stderr, "\n");
  }
  _Kdr_reg();
  _SlowCa_reg();
  _ampa_reg();
  _cad2_reg();
  _glutamate_reg();
  _h_reg();
  _h2_reg();
  _hh3_reg();
  _ih_reg();
  _it2_reg();
  _kap_reg();
  _kca_reg();
  _kdf_reg();
  _kdr2_reg();
  _km_reg();
}
