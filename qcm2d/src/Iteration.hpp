#include "Outputfile.h"
#include <Eigen/Sparse>

typedef Eigen::Triplet<double> T;
class Iteration {
	public:
			static void pu_pt(vvmArray &);
			static void pw_pt(vvmArray &);
			static void pth_pt(vvmArray &);
			static void ppi_pt(vvmArray &);
			static void pqv_pt(vvmArray &);
			static void pqc_pt(vvmArray &);
			static void pqr_pt(vvmArray &);
			static void condensation(vvmArray &, int, int); 	// condensation of qc by qv
			static void autoconversion(vvmArray &, int, int); 	// autoconversion of qc to qr
			static void accretion(vvmArray &, int, int); 		// accretion of qc by qr
			static void evaporation(vvmArray &, int, int); 	// evaporation of rain water
			static void heatflux(vvmArray &, int, int, int); 		// heat flux at surface 
			static void LeapFrog(vvmArray &);
};