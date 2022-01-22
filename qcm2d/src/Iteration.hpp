#include "Outputfile.h"
#include <Eigen/Sparse>

typedef Eigen::Triplet<double> T;
class Iteration {
	public:
			static void pu_pt(Array &);
			static void pw_pt(Array &);
			static void pth_pt(Array &);
			static void ppi_pt(Array &);
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