#include "Outputfile.hpp"
#include <Eigen/Sparse>

typedef Eigen::Triplet<double> T;
class Iteration {
	public:
			static void pzeta_pt(Array &);
			static void pth_pt(Array &);
			static void cal_w(Array &);
			static void cal_u(Array &);
			static void pqv_pt(Array &);
			static void pqc_pt(Array &);
			static void pqr_pt(Array &);
			static void condensation(Array &, int, int); 	// condensation of qc by qv
			static void autoconversion(Array &, int, int); 	// autoconversion of qc to qr
			static void accretion(Array &, int, int); 		// accretion of qc by qr
			static void evaporation(Array &, int, int); 	// evaporation of rain water
			static void heatflux(Array &, int, int, int); 		// heat flux at surface 
			static void LeapFrog(Array &);
};