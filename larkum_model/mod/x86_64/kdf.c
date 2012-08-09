/* Created by Language version: 6.2.0 */
/* NOT VECTORIZED */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib.h"
#undef PI
 
#include "md1redef.h"
#include "section.h"
#include "md2redef.h"

#if METHOD3
extern int _method3;
#endif

#undef exp
#define exp hoc_Exp
extern double hoc_Exp();
 
#define _threadargscomma_ /**/
#define _threadargs_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 static double *_p; static Datum *_ppvar;
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define gkbar _p[0]
#define vshift _p[1]
#define gk _p[2]
#define ninf _p[3]
#define ntau _p[4]
#define htau _p[5]
#define hinf _p[6]
#define ikd _p[7]
#define n _p[8]
#define h _p[9]
#define ek _p[10]
#define ik _p[11]
#define Dn _p[12]
#define Dh _p[13]
#define _g _p[14]
#define _ion_ek	*_ppvar[0]._pval
#define _ion_ik	*_ppvar[1]._pval
#define _ion_dikdv	*_ppvar[2]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 static int hoc_nrnpointerindex =  -1;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static int _hoc_rates();
 static int _hoc_trates();
 static int _mechtype;
extern int nrn_get_mechtype();
 static _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range();
 _prop = hoc_getdata_range(_mechtype);
 _p = _prop->param; _ppvar = _prop->dparam;
 ret(1.);
}
 /* connect user functions to hoc names */
 static IntFunc hoc_intfunc[] = {
 "setdata_kdf", _hoc_setdata,
 "rates_kdf", _hoc_rates,
 "trates_kdf", _hoc_trates,
 0, 0
};
 /* declare global and static user variables */
#define N N_kdf
 double N = 3;
#define hinf_b hinf_b_kdf
 double hinf_b = -11;
#define hinf_a hinf_a_kdf
 double hinf_a = -54;
#define htau_C6 htau_C6_kdf
 double htau_C6 = -55;
#define htau_C5 htau_C5_kdf
 double htau_C5 = 22;
#define htau_C4 htau_C4_kdf
 double htau_C4 = 48;
#define htau_C3 htau_C3_kdf
 double htau_C3 = -75;
#define htau_C2 htau_C2_kdf
 double htau_C2 = 810;
#define htau_C1 htau_C1_kdf
 double htau_C1 = 360;
#define ninf_b ninf_b_kdf
 double ninf_b = 14;
#define ninf_a ninf_a_kdf
 double ninf_a = -24;
#define ntau_C4 ntau_C4_kdf
 double ntau_C4 = 13;
#define ntau_C3 ntau_C3_kdf
 double ntau_C3 = -0.026;
#define ntau_C2 ntau_C2_kdf
 double ntau_C2 = 171;
#define ntau_C1 ntau_C1_kdf
 double ntau_C1 = 1.25;
#define ntau_slow2 ntau_slow2_kdf
 double ntau_slow2 = 0;
#define ntau_slow1 ntau_slow1_kdf
 double ntau_slow1 = 1;
#define q10 q10_kdf
 double q10 = 2.3;
#define tadj tadj_kdf
 double tadj = 0;
#define temp temp_kdf
 double temp = 23;
#define usetable usetable_kdf
 double usetable = 1;
#define vmax vmax_kdf
 double vmax = 100;
#define vmin vmin_kdf
 double vmin = -150;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "usetable_kdf", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "temp_kdf", "degC",
 "vmin_kdf", "mV",
 "vmax_kdf", "mV",
 "ntau_C1_kdf", "ms",
 "ntau_C3_kdf", "1/mV",
 "ntau_C4_kdf", "ms",
 "ninf_b_kdf", "mV",
 "htau_C1_kdf", "ms",
 "htau_C3_kdf", "mV",
 "htau_C4_kdf", "mV",
 "htau_C6_kdf", "mV",
 "hinf_a_kdf", "mV",
 "hinf_b_kdf", "mV",
 "gkbar_kdf", "pS/um2",
 "vshift_kdf", "mV",
 "gk_kdf", "pS/um2",
 "ntau_kdf", "ms",
 "htau_kdf", "ms",
 "ikd_kdf", "mA/cm2",
 0,0
};
 static double delta_t = 1;
 static double h0 = 0;
 static double n0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "temp_kdf", &temp_kdf,
 "q10_kdf", &q10_kdf,
 "vmin_kdf", &vmin_kdf,
 "vmax_kdf", &vmax_kdf,
 "ntau_C1_kdf", &ntau_C1_kdf,
 "ntau_C2_kdf", &ntau_C2_kdf,
 "ntau_C3_kdf", &ntau_C3_kdf,
 "ntau_C4_kdf", &ntau_C4_kdf,
 "ninf_a_kdf", &ninf_a_kdf,
 "ninf_b_kdf", &ninf_b_kdf,
 "htau_C1_kdf", &htau_C1_kdf,
 "htau_C2_kdf", &htau_C2_kdf,
 "htau_C3_kdf", &htau_C3_kdf,
 "htau_C4_kdf", &htau_C4_kdf,
 "htau_C5_kdf", &htau_C5_kdf,
 "htau_C6_kdf", &htau_C6_kdf,
 "hinf_a_kdf", &hinf_a_kdf,
 "hinf_b_kdf", &hinf_b_kdf,
 "ntau_slow1_kdf", &ntau_slow1_kdf,
 "ntau_slow2_kdf", &ntau_slow2_kdf,
 "N_kdf", &N_kdf,
 "tadj_kdf", &tadj_kdf,
 "usetable_kdf", &usetable_kdf,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(), nrn_init(), nrn_state();
 static void nrn_cur(), nrn_jacob();
 
static int _ode_count(), _ode_map(), _ode_spec(), _ode_matsol();
 
#define _cvode_ieq _ppvar[3]._i
 /* connect range variables in _p that hoc is supposed to know about */
 static char *_mechanism[] = {
 "6.2.0",
"kdf",
 "gkbar_kdf",
 "vshift_kdf",
 0,
 "gk_kdf",
 "ninf_kdf",
 "ntau_kdf",
 "htau_kdf",
 "hinf_kdf",
 "ikd_kdf",
 0,
 "n_kdf",
 "h_kdf",
 0,
 0};
 static Symbol* _k_sym;
 
static void nrn_alloc(_prop)
	Prop *_prop;
{
	Prop *prop_ion, *need_memb();
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 15, _prop);
 	/*initialize range parameters*/
 	gkbar = 2.5;
 	vshift = 0;
 	_prop->param = _p;
 	_prop->param_size = 15;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 4, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_k_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ek */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ik */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dikdv */
 
}
 static _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _update_ion_pointer(Datum*);
 _kdf_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("k", -10000.);
 	_k_sym = hoc_lookup("k_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
  hoc_register_dparam_size(_mechtype, 4);
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 kdf /home/torbjone/bedlewo/project/larkum_model/mod/x86_64/kdf.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double *_t_ninf;
 static double *_t_ntau;
 static double *_t_hinf;
 static double *_t_htau;
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static _modl_cleanup(){ _match_recurse=1;}
static _f_trates();
static rates();
static trates();
 
static int _ode_spec1(), _ode_matsol1();
 static _n_trates();
 static int _slist1[2], _dlist1[2];
 static int states();
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   trates ( _threadargscomma_ v - vshift ) ;
   ntau = ntau * ntau_slow1 + ntau_slow2 ;
   Dn = ( 1.0 - exp ( ( - dt * tadj ) / ( ntau ) ) ) * ( ninf - n ) / dt ;
   Dh = ( 1.0 - exp ( ( - dt * tadj ) / ( htau ) ) ) * ( hinf - h ) / dt ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 trates ( _threadargscomma_ v - vshift ) ;
 ntau = ntau * ntau_slow1 + ntau_slow2 ;
 Dn = Dn  / (1. - dt*( ( (( 1.0 - exp ( ( - dt * tadj ) / ( ntau ) ) ))*(( ( - 1.0 ) )) ) / dt )) ;
 Dh = Dh  / (1. - dt*( ( (( 1.0 - exp ( ( - dt * tadj ) / ( htau ) ) ))*(( ( - 1.0 ) )) ) / dt )) ;
}
 /*END CVODE*/
 static int states () {_reset=0;
 {
   trates ( _threadargscomma_ v - vshift ) ;
   ntau = ntau * ntau_slow1 + ntau_slow2 ;
    n = n + (1. - exp(dt*(( (( 1.0 - exp ( ( - dt * tadj ) / ( ntau ) ) ))*(( ( - 1.0 ) )) ) / dt)))*(- ( ( (( 1.0 - exp ( ( - dt * tadj ) / ( ntau ) ) ))*(( ninf )) ) / dt ) / ( ( (( 1.0 - exp ( ( - dt * tadj ) / ( ntau ) ) ))*(( ( - 1.0) )) ) / dt ) - n) ;
    h = h + (1. - exp(dt*(( (( 1.0 - exp ( ( - dt * tadj ) / ( htau ) ) ))*(( ( - 1.0 ) )) ) / dt)))*(- ( ( (( 1.0 - exp ( ( - dt * tadj ) / ( htau ) ) ))*(( hinf )) ) / dt ) / ( ( (( 1.0 - exp ( ( - dt * tadj ) / ( htau ) ) ))*(( ( - 1.0) )) ) / dt ) - h) ;
   }
  return 0;
}
 static double _mfac_trates, _tmin_trates;
 static _check_trates();
 static _check_trates() {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  static double _sav_dt;
  static double _sav_celsius;
  static double _sav_temp;
  static double _sav_ntau_C1;
  static double _sav_ntau_C2;
  static double _sav_ntau_C3;
  static double _sav_ntau_C4;
  static double _sav_ninf_a;
  static double _sav_ninf_b;
  static double _sav_htau_C1;
  static double _sav_htau_C2;
  static double _sav_htau_C3;
  static double _sav_htau_C4;
  static double _sav_htau_C5;
  static double _sav_htau_C6;
  static double _sav_hinf_a;
  static double _sav_hinf_b;
  if (!usetable) {return;}
  if (_sav_dt != dt) { _maktable = 1;}
  if (_sav_celsius != celsius) { _maktable = 1;}
  if (_sav_temp != temp) { _maktable = 1;}
  if (_sav_ntau_C1 != ntau_C1) { _maktable = 1;}
  if (_sav_ntau_C2 != ntau_C2) { _maktable = 1;}
  if (_sav_ntau_C3 != ntau_C3) { _maktable = 1;}
  if (_sav_ntau_C4 != ntau_C4) { _maktable = 1;}
  if (_sav_ninf_a != ninf_a) { _maktable = 1;}
  if (_sav_ninf_b != ninf_b) { _maktable = 1;}
  if (_sav_htau_C1 != htau_C1) { _maktable = 1;}
  if (_sav_htau_C2 != htau_C2) { _maktable = 1;}
  if (_sav_htau_C3 != htau_C3) { _maktable = 1;}
  if (_sav_htau_C4 != htau_C4) { _maktable = 1;}
  if (_sav_htau_C5 != htau_C5) { _maktable = 1;}
  if (_sav_htau_C6 != htau_C6) { _maktable = 1;}
  if (_sav_hinf_a != hinf_a) { _maktable = 1;}
  if (_sav_hinf_b != hinf_b) { _maktable = 1;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_trates =  vmin ;
   _tmax =  vmax ;
   _dx = (_tmax - _tmin_trates)/199.; _mfac_trates = 1./_dx;
   for (_i=0, _x=_tmin_trates; _i < 200; _x += _dx, _i++) {
    _f_trates(_x);
    _t_ninf[_i] = ninf;
    _t_ntau[_i] = ntau;
    _t_hinf[_i] = hinf;
    _t_htau[_i] = htau;
   }
   _sav_dt = dt;
   _sav_celsius = celsius;
   _sav_temp = temp;
   _sav_ntau_C1 = ntau_C1;
   _sav_ntau_C2 = ntau_C2;
   _sav_ntau_C3 = ntau_C3;
   _sav_ntau_C4 = ntau_C4;
   _sav_ninf_a = ninf_a;
   _sav_ninf_b = ninf_b;
   _sav_htau_C1 = htau_C1;
   _sav_htau_C2 = htau_C2;
   _sav_htau_C3 = htau_C3;
   _sav_htau_C4 = htau_C4;
   _sav_htau_C5 = htau_C5;
   _sav_htau_C6 = htau_C6;
   _sav_hinf_a = hinf_a;
   _sav_hinf_b = hinf_b;
  }
 }

 static trates(double _lv){ _check_trates();
 _n_trates(_lv);
 return;
 }

 static _n_trates(double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 _f_trates(_lv); return; 
}
 _xi = _mfac_trates * (_lv - _tmin_trates);
 _i = (int) _xi;
 if (_xi <= 0.) {
 ninf = _t_ninf[0];
 ntau = _t_ntau[0];
 hinf = _t_hinf[0];
 htau = _t_htau[0];
 return; }
 if (_i >= 199) {
 ninf = _t_ninf[199];
 ntau = _t_ntau[199];
 hinf = _t_hinf[199];
 htau = _t_htau[199];
 return; }
 _theta = _xi - (double)_i;
 ninf = _t_ninf[_i] + _theta*(_t_ninf[_i+1] - _t_ninf[_i]);
 ntau = _t_ntau[_i] + _theta*(_t_ntau[_i+1] - _t_ntau[_i]);
 hinf = _t_hinf[_i] + _theta*(_t_hinf[_i+1] - _t_hinf[_i]);
 htau = _t_htau[_i] + _theta*(_t_htau[_i+1] - _t_htau[_i]);
 }

 
static int  _f_trates (  _lv )  
	double _lv ;
 {
   rates ( _threadargscomma_ _lv ) ;
    return 0; }
 
static int _hoc_trates() {
  double _r;
    _r = 1.;
 trates (  *getarg(1) ) ;
 ret(_r);
}
 
static int  rates (  _lv )  
	double _lv ;
 {
   double _lntau_C2_50 , _lntau_C3_50 ;
 _lntau_C2_50 = ntau_C2 ;
   _lntau_C3_50 = ntau_C3 ;
   if ( _lv > - 50.0 ) {
     _lntau_C2_50 = ntau_C4 ;
     _lntau_C3_50 = - ntau_C3 ;
     }
   ntau = ntau_C1 + _lntau_C2_50 * exp ( - _lv * _lntau_C3_50 ) ;
   ntau = ntau ;
   htau = htau_C1 + ( htau_C2 + htau_C5 * ( _lv - htau_C6 ) ) * exp ( - pow( ( ( _lv - htau_C3 ) / htau_C4 ) , 2.0 ) ) ;
   ninf = 1.0 / ( 1.0 + exp ( ( ninf_a - _lv ) / ninf_b ) ) ;
   hinf = 1.0 / ( 1.0 + exp ( ( hinf_a - _lv ) / hinf_b ) ) ;
    return 0; }
 
static int _hoc_rates() {
  double _r;
   _r = 1.;
 rates (  *getarg(1) ) ;
 ret(_r);
}
 
static int _ode_count(_type) int _type;{ return 2;}
 
static int _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ek = _ion_ek;
     _ode_spec1 ();
  }}
 
static int _ode_map(_ieq, _pv, _pvdot, _pp, _ppd, _atol, _type) int _ieq, _type; double** _pv, **_pvdot, *_pp, *_atol; Datum* _ppd; { 
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 2; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static int _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ek = _ion_ek;
 _ode_matsol1 ();
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_k_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_k_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_k_sym, _ppvar, 2, 4);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  h = h0;
  n = n0;
 {
   trates ( _threadargscomma_ v - vshift ) ;
   n = ninf ;
   h = hinf ;
   tadj = pow( q10 , ( ( celsius - temp ) / 10.0 ) ) ;
   }
  _sav_indep = t; t = _save;

}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
  ek = _ion_ek;
 initmodel();
 }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   gk = gkbar * h * ( tadj ) * pow( n , N ) ;
   ik = ( 1e-4 ) * gk * ( v - ek ) ;
   ikd = ik ;
   }
 _current += ik;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
  ek = _ion_ek;
 _g = _nrn_current(_v + .001);
 	{ double _dik;
  _dik = ik;
 _rhs = _nrn_current(_v);
  _ion_dikdv += (_dik - ik)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ik += ik ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type){
 double _break, _save;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 _break = t + .5*dt; _save = t;
 v=_v;
{
  ek = _ion_ek;
 { {
 for (; t < _break; t += dt) {
 error =  states();
 if(error){fprintf(stderr,"at line 95 in file kdf.mod:\n      SOLVE states METHOD cnexp\n"); nrn_complain(_p); abort_run(error);}
 
}}
 t = _save;
 } }}

}

static terminal(){}

static _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(n) - _p;  _dlist1[0] = &(Dn) - _p;
 _slist1[1] = &(h) - _p;  _dlist1[1] = &(Dh) - _p;
   _t_ninf = makevector(200*sizeof(double));
   _t_ntau = makevector(200*sizeof(double));
   _t_hinf = makevector(200*sizeof(double));
   _t_htau = makevector(200*sizeof(double));
_first = 0;
}
