
#include "SiStripGeomFTD.h"

// Include standard header files
#include "Colours.h"
#include "PhysicalConstants.h"

#include <cstdlib>
#include <iomanip>
#include <algorithm>
#include <cmath>

// Include CLHEP header files
#include <CLHEP/Vector/EulerAngles.h>
#include <CLHEP/Vector/Rotation.h>

// Include Gear header files
#include <gear/BField.h>
//#include <gear/GearParameters.h>
#include <gearimpl/FTDParametersImpl.h>
#include <gearimpl/FTDLayerLayoutImpl.h>
#include <gearimpl/Vector3D.h>

// Include Marlin
#include <marlin/Global.h>
#include <streamlog/streamlog.h>

// Include Encoder stuff
#include <UTIL/ILDConf.h>



// FIXME: Maybe to be changed in gear xml. By the moment
//        the xml gives mm as distance unit whereas this
//        module uses cm.

//template <class T> void convertunits(std::vector<T> & v2convert);

//FIXME: DEPRECATED
template <class T>
void convertunits(std::vector<T> & vectoconvert, const T unit)
{
	for(unsigned int i=0; i < vectoconvert.size(); i++)
	{
		vectoconvert[i] = vectoconvert[i]*unit;
	}
}


namespace sistrip 
{
SiStripGeomFTD::SiStripGeomFTD(const std::string & detector):
	SiStripGeom(detector),_ftdParams(0),_ftdLayer(0)
{
	if( _gearType != "FTD")
	{
		streamlog_out(ERROR) << "Wrong instantiation of the class. " <<
			"See SiStripBuilder class, should have an incoherence!" <<
			std::endl;
		exit(0);
	}
}

SiStripGeomFTD::~SiStripGeomFTD()
{
}


// Method initializing this class - reads Gear parameters from XML file
void SiStripGeomFTD::initGearParams()
{

	// BField: needed for every tracker detector
	try 
	{
		const gear::BField & bField = marlin::Global::GEAR->getBField();
		
		gear::Vector3D auxvec(0.0,0.0,0.0);
		
		_magField.setX( (bField.at( auxvec )).x() * T);
		_magField.setY( (bField.at( auxvec )).y() * T);
		_magField.setZ( (bField.at( auxvec )).z() * T);
	}
	catch(gear::UnknownParameterException& e) 
	{
		std::cout << "No magnetic field found in gear file!" << std::endl;
	}

	//------Get the geometry from the gear file-----//
	try
	{
		_ftdParams = const_cast<gear::FTDParameters *>( &(marlin::Global::GEAR->getFTDParameters()) );
	}
	catch(gear::UnknownParameterException& e)
	{
		std::cout << "No FTD found in gear file" << std::endl;
		exit(-1);
	}

	try
	{
		_ftdLayer = const_cast<gear::FTDLayerLayout *>( &_ftdParams->getFTDLayerLayout() );
	}
	catch(gear::UnknownParameterException  & e)
	{
		std::cout << "Unexpected Error! No FTDLayerLayout found in FTDParameters class" << std::endl;
		exit(-1);
	}

	//FIXME: Data members of the SiStripGeom, this is needed ?? If SiStripGeom is the interface.. yes
	_numberOfLayers  = _ftdLayer->getNLayers();
	// Layers goes from 0,..., 2N-1
	_layerType.reserve(2*_numberOfLayers);
	_layerZ.reserve(2*_numberOfLayers);
	_layerRadius.reserve(2*_numberOfLayers); 
	_layerPhi0.reserve(2*_numberOfLayers);
	_layerHalfPhi.reserve(2*_numberOfLayers);
	_layerPetalOpAngle.reserve(2*_numberOfLayers);
	_layerOuterRadius.reserve(2*_numberOfLayers);
	_sensorThick.reserve(2*_numberOfLayers);
	_sensorWidth.reserve(2*_numberOfLayers);
	_sensorWidth2.reserve(2*_numberOfLayers);
	_sensorLength.reserve(2*_numberOfLayers);
	_layerRealID.reserve(2*_numberOfLayers);
	_ladderOffsetZ.reserve(2*_numberOfLayers);
	_ladderZOffsetSign0.reserve(2*_numberOfLayers);
	_numberOfLadders.reserve(2*_numberOfLayers);
	_ladderLength.reserve(2*_numberOfLayers);
	for(int i = 0; i < _numberOfLayers; i++)
	{
		_layerType.push_back(_ftdLayer->getSensorType(i));
		_layerZ.push_back(_ftdLayer->getZposition(i)*mm);
		_layerRadius.push_back( _ftdLayer->getSensitiveRinner(i)*mm);
		_layerPhi0.push_back(_ftdLayer->getPhi0(i));
		_layerHalfPhi.push_back(_ftdLayer->getPhiHalfDistance(i));
		_layerPetalOpAngle.push_back(_ftdLayer->getPhiHalfDistance(i));
		_layerOuterRadius.push_back(_ftdLayer->getSensitiveRinner(i)+
				_sensorLength[i]*mm);
		_layerRealID.push_back(i+1);      

		_sensorThick.push_back(_ftdLayer->getSensitiveThickness(i)*mm);
		_sensorWidth.push_back(_ftdLayer->getSensitiveLengthMax(i)*mm);// x-direction
		_sensorWidth2.push_back(_ftdLayer->getSensitiveLengthMin(i)*mm);// x-direction
		_sensorLength.push_back(_ftdLayer->getSensitiveWidth(i)*mm);   // y-direction
		_ladderZOffsetSign0.push_back(_ftdLayer->getZoffsetSign0(i));
		_ladderOffsetZ.push_back(_ftdLayer->getZoffset(i)*mm);
		_numberOfLadders.push_back(_ftdLayer->getNPetals(i));
		_ladderLength.push_back(_ftdLayer->getSupportWidth(i)*mm);
	}
	// Negative layers
	const unsigned int ridsize = _layerRealID.size();
	for(unsigned int i=0; i < ridsize; i++)
	{
		_layerRealID.push_back(-1*_layerRealID.at(i));
	}
	_layerType.insert(_layerType.end(),_layerType.begin(),_layerType.end());
	_layerZ.insert(_layerZ.end(),_layerZ.begin(),_layerZ.end());
	_layerRadius.insert(_layerRadius.end(),_layerRadius.begin(),_layerRadius.end());
	_layerPhi0.insert(_layerPhi0.end(),_layerPhi0.begin(),_layerPhi0.end());
	_layerHalfPhi.insert(_layerHalfPhi.end(),_layerHalfPhi.begin(),_layerHalfPhi.end());
	_layerPetalOpAngle.insert(_layerPetalOpAngle.end(),_layerPetalOpAngle.begin(),
			_layerPetalOpAngle.end());
	_layerOuterRadius.insert(_layerOuterRadius.end(),_layerOuterRadius.begin(),
			_layerOuterRadius.end());
	_sensorThick.insert(_sensorThick.end(),_sensorThick.begin(),_sensorThick.end());
	_sensorWidth.insert(_sensorWidth.end(),_sensorWidth.begin(),_sensorWidth.end());
	_sensorWidth2.insert(_sensorWidth2.end(),_sensorWidth2.begin(),_sensorWidth2.end());
	_sensorLength.insert(_sensorLength.end(),_sensorLength.begin(),
			_sensorLength.end());
	_ladderOffsetZ.insert(_ladderOffsetZ.end(),_ladderOffsetZ.begin(),
			_ladderOffsetZ.end());
	_ladderZOffsetSign0.insert(_ladderZOffsetSign0.end(),_ladderZOffsetSign0.begin(),
			_ladderZOffsetSign0.end());
	_numberOfLadders.insert(_numberOfLadders.end(),_numberOfLadders.begin(),
			_numberOfLadders.end());
	_ladderLength.insert(_ladderLength.end(),_ladderLength.begin(),
			_ladderLength.end());
	
	// FIXME: Not for the gear--> must be a processor input parameter
	//_sensorPitchInFront  = std::vector<double>(2*_numberOfLayers,50.0 * um); //FIXME: HARDCODED
	_sensorPitchInFront  = std::vector<double>(2*_numberOfLayers,50.0 * um); //FIXME: HARDCODED
	_sensorNStripsInFront.reserve(2*_numberOfLayers);
	_sensorPitchInRear     = std::vector<double>(2*_numberOfLayers,50.0*um); // FIXME
	_sensorNStripsInRear.reserve(2*_numberOfLayers);
	for(int i = 0; i < _numberOfLayers; i++)
	{
		// 50 um in the middle
		const double xmaxsensorup = _ftdLayer->getSensitiveLengthMax(i)*mm;
		const double ymiddlesensorup  = _ftdLayer->getSensitiveWidth(i)*mm/2.0;
		const double widthmiddle = xmaxsensorup-
			(2.0*tan(_ftdLayer->getPhiHalfDistance(i))*ymiddlesensorup);
		
		_sensorNStripsInFront.push_back( (int)(widthmiddle/_sensorPitchInFront[i])+1 );
		_sensorNStripsInRear.push_back( (int)(widthmiddle/_sensorPitchInRear[i])+1 );
	}
	_sensorNStripsInFront.insert(_sensorNStripsInFront.end(),_sensorNStripsInFront.begin(),
			_sensorNStripsInFront.end());
	_sensorNStripsInRear.insert(_sensorNStripsInRear.end(),_sensorNStripsInRear.begin(),
		_sensorNStripsInRear.end());

	/*_sensorPitchInRear     = std::vector<double>(2*_numberOfLayers,50.0*um); // FIXME
	_sensorNStripsInRear.reserve(2*_numberOfLayers);
	for(int i = 0; i < _numberOfLayers; i++)
	{
		const double widthsensor = _ftdLayer->getSensitiveWidth(i)*mm;
		_sensorNStripsInRear.push_back( (int)(widthsensor/_sensorPitchInRear[i])+1 );
	}
	_sensorNStripsInRear.insert(_sensorNStripsInRear.end(),_sensorNStripsInRear.begin(),
		_sensorNStripsInRear.end());*/
}


void SiStripGeomFTD::updateCanonicalCellID(const int & cellID, const int & stripType,
		const int & stripID, UTIL::BitField64 * cellEnc)
{
	(*cellEnc)["subdet"]=ILDDetID::FTD;
	short int layer=-1;
	short int module=-1;
	short int sensor=-1;

	decodeCellID(layer,module,sensor,cellID);

	int realLayer = getLayerRealID(layer);

	(*cellEnc)["side"]=abs(realLayer)/realLayer;
	(*cellEnc)["layer"]=abs(realLayer);
	(*cellEnc)["module"]=module;
	(*cellEnc)["sensor"]=sensor;
	(*cellEnc)["stripType"]=stripType;
	(*cellEnc)["stripID"]=stripID;	
}


// FIXME: TO BE DEPRECATED!!
std::map<std::string,short int> SiStripGeomFTD::cellIDDecProv(EVENT::SimTrackerHit * & simHit)
{
	// Encoded disks: 0,...,6  positives
	//                7,...,13 negatives
	const static int SHIFT_LAYER =0;   // DISKS
	const static int SHIFT_LADDER=9;   // LADDER
	const static int SHIFT_SENSOR=17;  // SENSOR
	
	const static unsigned int MASK_LAYER =(unsigned int)0x000001FF; 
 	const static unsigned int MASK_LADDER=(unsigned int)0x0001FE00;
	const static unsigned int MASK_SENSOR=(unsigned int)0x01FE0000;
	const static unsigned int NBITS_LAYER=(unsigned int)0x000001FF;


	const unsigned int code = simHit->getCellID();

        unsigned int layerid = (code & MASK_LAYER) >> SHIFT_LAYER;
        unsigned int ladderid= (code & MASK_LADDER) >> SHIFT_LADDER;
        unsigned int sensorid= (code & MASK_SENSOR) >> SHIFT_SENSOR;

        // Checking the negative Z
        if( layerid > 7 )
        {
                // 
                layerid = (~layerid & NBITS_LAYER)+1;
		layerid += 7;
        }
	std::map<std::string,short int> id;
        id["layer"] = (layerid-1);
        id["ladder"]= ladderid-1;
        id["sensor"]= sensorid;
	
	return id;
}

// Added
double SiStripGeomFTD::getLadderOffsetX(const short int & layerID) const
{

	if( _sensorWidth.size() < (unsigned int)layerID )
	{
		streamlog_out(ERROR) << 
			"SiStripGeomFTD::getLadderOffsetX - layerID: " << 
			layerID << " out of range!!!" << std::endl;
	}
	
	return _sensorWidth[layerID]/2.0;
}

// TRANSFORMATION METHODS - GLOBAL REF. SYSTEM

// The Local Reference Frame system (LRF) is defined in the way that:
//  - is defined positive
//  - the electric drift field is oriented in the negative direction
//    of the X-axis
//  - all the hits inside a sensors have positive local coordinates
//  - the strips collecting holes (Single Side as Double side sensors)
//    are faced to the IP
//  - the strips collecting holes are in the sensor plane X=0 (the holes
//    are collected by the strips oriented in Z-direction).
// These requerimemts make placing the LRF:
//
//                             --------
//                             |      |              Face near IP
//    1 and 2 sensors           |    |                
//                               |__| X <--- 
//
//
//                             --------
//                             |      |
//    3 and 4 sensors           |    |               Face far IP
//                        ---> X |__|  
//


//
// Method transforming given point from global ref. system to local ref. system
// (parameters: diskID, petalID, sensorID and space point in global ref. system)
//
CLHEP::Hep3Vector SiStripGeomFTD::transformPointToLocal(short int diskID, short int petalID, short int sensorID, const CLHEP::Hep3Vector & globalPoint)
{
	// Initialize local point
	CLHEP::Hep3Vector localPoint(globalPoint);

	// Phi Angle of the Petal (on the centroid of the petal with respect the x-axis)
	const double phi0 = _ftdLayer->getPhiPetalCd(diskID,petalID);
	
	// 
	// RECALL: ftd gear have the convention mm (distance), here we have cm
	// 

	// Extract the Z of the sensor: diferent position depending the sensor
	const double zsensorCd = _ftdLayer->getSensitiveZposition(diskID,petalID,sensorID)*mm;
	const int zsign = (int)(zsensorCd/fabs(zsensorCd));
	const double sensorthickness = _ftdLayer->getSensitiveThickness(diskID)*mm;
	// Sensor 3 and 4: Displacing to the trapezoid farest the IP
	double zsensor = zsensorCd+zsign*sensorthickness/2.0;
	// Sensor 1 and 2: Displacing to the trapezoid facing the IP
	if( sensorID < 3 )
	{
		zsensor = zsensorCd-zsign*sensorthickness/2.0;
	}
	
	// And the X and Y CentroiD position: over the smallest side of the trapezoid
	const int xsign = (int)(fabs(globalPoint.getX())/globalPoint.getX());
	const double xsensorCd = xsign*(_ftdLayer->getSensitiveRinner(diskID)*mm)*
		fabs(cos(phi0));

	const int ysign = (int)(fabs(globalPoint.getY())/globalPoint.getY());
	const double ysensorCd = ysign*(_ftdLayer->getSensitiveRinner(diskID)*mm)*
		fabs(sin(phi0));

	// The (0,0,0) position of LFR
	CLHEP::Hep3Vector localOrigin(xsensorCd, ysensorCd, zsensor);

	// Translating the globalPoint to the local Origin
	localPoint -= localOrigin;
	
	// Perform rotation to get the system local: 
	//  Z-positives sensors 1,2 |
	//  Z-negatives sensors 3,4 |-- Same transformation
	//
	//  Z-positives sensors 3,4 | ---> La antigua
	//  Z-negatives sensors 1,2 |-- Same transformation
	double rotZangle = -phi0;
	double rotYangle = -M_PI/2.0;
	if( (sensorID < 3 && zsign > 0) || 
			(sensorID > 2 && zsign < 0) )
	{
		rotZangle = M_PI-phi0;
		rotYangle = M_PI/2.0;
	}
	localPoint.rotateZ(rotZangle);
	localPoint.rotateY(rotYangle);

	// Avoiding X,Y and Z negative values--> displacing from the
        // centroid to the edge
        const double longtrapezoidedge = _ftdLayer->getSensitiveLengthMax(diskID)*mm;
	localPoint += CLHEP::Hep3Vector(0.0,longtrapezoidedge/2.,0.0); 
	
	streamlog_out(DEBUG4) << 
		"============================================\n" 
		<< "SiStripGeomFTD::transformPointToLocal \n" 
		<< " Sensor ID (C-vector style): \n "
		<< "   DISK: "<<diskID<<" PETAL: "<<petalID<<" SENSOR: " << sensorID << "\n"
		<< "   Petal Phi: " << phi0*180.0/M_PI << "\n"
		<< std::setprecision(3) 
		<< " Origen of Local frame [mm]: (" << xsensorCd/mm <<","
			<< ysensorCd/mm << "," << zsensor/mm << ")\n"
		<< " Hit (Global ref. frame) [mm]:" << globalPoint/mm << "\n"
		<< " Hit (Local ref. frame)  [mm]:" << localPoint/mm << "\n" 
		<< " Maximum dimensions: x < " << sensorthickness/mm 
			<<  ", y < " << longtrapezoidedge/mm << ", z < " 
			<<  _ftdLayer->getSensitiveWidth(diskID)*mm/mm << " [mm]\n"
		<< "============================================" << std::endl;
	
	// Return space point in local ref. system
	return localPoint;
}

//
// Method transforming given vector from global ref. system to local ref. system
// (parameters: layerID, ladderID, sensorID and vector in global ref. system)
//
CLHEP::Hep3Vector SiStripGeomFTD::transformVecToLocal(short int diskID, short int petalID,
		short int sensorID, const CLHEP::Hep3Vector & globalVec)
{
	// Initialize local vector
	CLHEP::Hep3Vector localVec(globalVec);

	const double phi0 = _ftdLayer->getPhiPetalCd(diskID,petalID);
	// Extract z-sign and 
	const int realLayerID = getLayerRealID(diskID);
	const int zsign = abs(realLayerID)/realLayerID;
	// Sensors 3,4 z-positive and 1,2 z-negative
	double rotZangle = -phi0;
	double rotYangle = -M_PI/2.0;
	// Sensors 1,2 z-positive and 3,4 z-negative
	if( (sensorID < 3 && zsign > 0) || 
			(sensorID > 2 && zsign < 0) )
	{
		rotZangle = M_PI-phi0;
		rotYangle = M_PI/2.0;
	}
	localVec.rotateZ(rotZangle);
	localVec.rotateY(rotYangle);
	
	// Return vector in local ref. system
	return localVec;
}

//
// Method transforming given matrix 3x3 from global ref. system to local ref. system
// (parameters: layerID, ladderID, sensorID and matrix in global ref. system)
//
CLHEP::HepMatrix SiStripGeomFTD::transformMatxToLocal(short int layerID, short int ladderID, const CLHEP::HepMatrix & globalMatrix)
{
	// FIXME: NOT DONE YET
	std::cout << "SiStripGeomFTD::transformMatxToLocal -- METHOD NOT IMPLEMENTED" 
		<< std::endl;
	exit(0);
	// Initialize local matrix 3x3 to zero values
	CLHEP::HepMatrix localMatrix(3,3,0);
	
	// Initialize rotation matrices: R, R^T (transposition)
	CLHEP::HepMatrix rotMatrix(3,3,0);
	CLHEP::HepMatrix rotMatrixT(3,3,0);

	// Gear type: VXD
	   if (_gearType == "VXD") {

	      // Calculate rotation angles
	      double theta = getLadderTheta(layerID); // ALWAYS M_PI/2
	      double phi   = getLadderPhi(layerID, ladderID);

	      // Calculate rotation matrices - help
	      CLHEP::HepRotation rotMatrixZ(CLHEP::Hep3Vector(0,0,1),-phi);
	      CLHEP::HepRotation rotMatrixY(CLHEP::Hep3Vector(0,1,0),+theta);
	      CLHEP::HepRotation rotMatrixHelp(rotMatrixY*rotMatrixZ);

	      CLHEP::HepRotation rotMatrixHelpT(rotMatrixHelp.inverse());

	      for (int i=0; i<3; i++) {
	         for (int j=0; j<3; j++) {
	            rotMatrix[i][j]  = rotMatrixHelp[i][j];
	            rotMatrixT[i][j] = rotMatrixHelpT[i][j];
	         }
	      }
	   }

	   // Gear type: unknown - error
	   else {
	      streamlog_out(ERROR) << "Unknown gear type!"
	                           << std::endl;

	      exit(0);
	   }

   // Transform given matrix - rotation wrt global ref. system
   localMatrix = rotMatrix*globalMatrix*rotMatrixT;

   // Return matrix in local ref. system
   return localMatrix;
}


// TRANSFORMATION METHODS - LOCAL REF. SYSTEM

//
// Method transforming given point from local ref. system to global ref. system
// (parameters: layerID, ladderID, sensorID and space point in local ref. system)
//
CLHEP::Hep3Vector SiStripGeomFTD::transformPointToGlobal(short int diskID, 
		short int petalID, short int sensorID, const CLHEP::Hep3Vector & localPoint)
{
	// Initialize global point
	CLHEP::Hep3Vector globalPoint(localPoint);

	// Calculate rotation angles
	double theta = getLadderTheta(diskID);  // ALWAYS M_PI/2.0
	double phi   = getLadderPhi(diskID, petalID);
	
	// Taking account of the z-side
	const int realLayerID = getLayerRealID(diskID);
	const int zsign = abs(realLayerID)/realLayerID;
	double rotZangle = -phi;
	double rotYangle = -theta;
	if( (sensorID < 3 && zsign > 0) || 
			(sensorID > 2 && zsign < 0) )
	{
		rotZangle = M_PI-phi;
		rotYangle = theta;
	}
	
	// Perform translation - to the center of a petal (in local frame)
	// Note that the center of the petal we decide defining in the back face
	// of the sensor
	const double xlocalCd = 0.0; // Already in the back face_ftdLayer->getSensitiveThickness(diskID)*mm;
	const double ylocalCd = _ftdLayer->getSensitiveLengthMax(diskID)/2.0*mm;
	const double zlocalCd = 0.0; // already in the low edge of the petal (see ysensorCd)
	
	globalPoint -= CLHEP::Hep3Vector(xlocalCd,ylocalCd,zlocalCd);

	// Perform rotation - with respect to the local centre of a ladder
	globalPoint.rotateY(-rotYangle);
	globalPoint.rotateZ(-rotZangle);
	// Now we have the vector in global coordinates from the center of the petal
	// to the hit 
	
	// Find (0,0,0) position of local coordinate system (in global coord.)
	// Extract the Z of the sensor
	const double zsensorCd = _ftdLayer->getSensitiveZposition(diskID,petalID,sensorID)*mm;
	const double sensorthickness = _ftdLayer->getSensitiveThickness(diskID)*mm;
	// Sensor 3 and 4: Displacing to the trapezoid farest the IP
	double zsensor = zsensorCd-zsign*sensorthickness/2.0;
	// Sensor 1 and 2: Displacing to the trapezoid facing the IP
	if( sensorID < 3 )
	{
		zsensor = zsensorCd+zsign*sensorthickness/2.0;
	}
	/* OLD
	// Sensors 3 and 4: Displacing to the backed the IP face
	double zsensorback = zsensor+zsign*sensorthickness/2.0;
	if( sensorID < 3 )
	{
		zsensorback = zsensor-zsign*sensorthickness/2.0;
	} END OLD*/
	// And the X and Y CentroiD position: over the smallest side of the trapezoid
	const double xsensorCd = _ftdLayer->getSensitiveRinner(diskID)*mm*cos(phi);
	const double ysensorCd = _ftdLayer->getSensitiveRinner(diskID)*mm*sin(phi);
	CLHEP::Hep3Vector localOrigin(xsensorCd,ysensorCd,zsensor);
	
	// Perform translation - to the global system
	globalPoint += localOrigin;
	
	streamlog_out(DEBUG4) << 
		"============================================\n" 
		<< "SiStripGeomFTD::transformPointToGlobal \n" 
		<< " Sensor ID (C-vector style): \n "
		<< "   DISK: "<<diskID<<" PETAL: "<<petalID<<" SENSOR: " << sensorID << "\n"
		<< "   Petal Phi: " << phi*180.0/M_PI << "\n"
		<< std::setprecision(3) 
		<< " Origen of Local frame [mm]: (" << xsensorCd/mm <<","
			<< ysensorCd/mm << "," << zsensor/mm << ")\n"
		<< " Hit (Global ref. frame) [mm]:" << globalPoint/mm << "\n"
		<< " Hit (Local ref. frame)  [mm]:" << localPoint/mm << "\n" 
		<< "============================================" << std::endl;
	
	// Return space point in global ref. system
	return globalPoint;
}

//
// Method transforming given vector from local ref. system to global ref. system
// (parameters: layerID, ladderID, sensorID and vector in local ref. system)
//
CLHEP::Hep3Vector SiStripGeomFTD::transformVecToGlobal(short int diskID, short int petalID, 
		short int sensorID, const CLHEP::Hep3Vector & localVec)
{
	// Initialize global vector
	CLHEP::Hep3Vector globalVec(localVec);

	
	// Calculate rotation angles
	double theta = getLadderTheta(diskID); // ALWAYS PI/2
	double phi   = getLadderPhi(diskID, petalID);
	
	const int realLayerID = getLayerRealID(diskID);
	const int zsign = abs(realLayerID)/realLayerID;
	double rotZangle = -phi;
	double rotYangle = -theta;
	if( (sensorID < 3 && zsign > 0) || 
			(sensorID > 2 && zsign < 0) )
	{
		rotZangle = M_PI-phi;
		rotYangle = theta;
	}
	globalVec.rotateY(-rotYangle);
	globalVec.rotateZ(-rotZangle);
	
	// Return vector in global ref. system
	return globalVec;
}

//
// Method transforming given matrix 3x3 from local ref. system to global ref. system
// (parameters: layerID, ladderID, sensorID and matrix in local ref. system)
//
CLHEP::HepMatrix SiStripGeomFTD::transformMatxToGlobal(short int layerID, short int ladderID, const CLHEP::HepMatrix & localMatrix)
{
	// FIXME: NOT DONE YET
	std::cout << "SiStripGeomFTD::transformMatxToGlobal -- METHOD NOT IMPLEMENTED" 
		<< std::endl;
	exit(0);
   // Initialize local matrix 3x3 to zero values
   CLHEP::HepMatrix globalMatrix(3,3,0);

   // Initialize rotation matrices: R, R^T (transposition)
   CLHEP::HepMatrix rotMatrix(3,3,0);
   CLHEP::HepMatrix rotMatrixT(3,3,0);

   // Gear type: VXD
   if (_gearType == "VXD") {

      // Calculate rotation angles
      double theta = getLadderTheta(layerID);
      double phi   = getLadderPhi(layerID, ladderID);

      // Calculate rotation matrices - help
      CLHEP::HepRotation rotMatrixY(CLHEP::Hep3Vector(0,1,0),-theta);
      CLHEP::HepRotation rotMatrixZ(CLHEP::Hep3Vector(0,0,1),+phi);
      CLHEP::HepRotation rotMatrixHelp(rotMatrixZ*rotMatrixY);

      CLHEP::HepRotation rotMatrixHelpT(rotMatrixHelp.inverse());

      for (int i=0; i<3; i++) {
         for (int j=0; j<3; j++) {
            rotMatrix[i][j]  = rotMatrixHelp[i][j];
            rotMatrixT[i][j] = rotMatrixHelpT[i][j];
         }
      }
   }

   // Gear type: unknown - error
   else {
      streamlog_out(ERROR) << "Unknown gear type!"
                           << std::endl;

      exit(0);
   }

   // Transform given matrix - rotation wrt local ref. system
   globalMatrix = rotMatrix*localMatrix*rotMatrixT;

   // Return matrix in global ref. system
   return globalMatrix;
}


// OTHER METHODS - GLOBAL REF. SYSTEM

//
// Get info whether the given point is inside of Si sensor (parameters: layerID,
// space point in local ref. system)
//
bool SiStripGeomFTD::isPointInsideSensor (short int diskID, short int petalID, 
		short int sensorID, const CLHEP::Hep3Vector & point) const
{
	// Angle defining the petals
	const double phi0 = _ftdLayer->getPhiHalfDistance(diskID);

	// Calculate the maximum of each axis
	const double xmax = _ftdLayer->getSensitiveThickness(diskID)*mm;
	const double ygap = (_ftdLayer->getSensitiveWidth(diskID)*mm - point.getZ())*tan(phi0);
	const double ymax = _ftdLayer->getSensitiveLengthMax(diskID)*mm-ygap;
	const double zmax = _ftdLayer->getSensitiveWidth(diskID)*mm;

	bool isIn = true;
	// Boundary set +- epsilon
        if( (point.getX() > xmax+EPS*um) || (point.getX() < -EPS*um) )
	{
		streamlog_out(ERROR) << " position out of sensor! X [um] = "
			<< point.getX()/um << " (Maximum x: " << xmax/um << ")" << std::endl;
		isIn = false;
	}

	if( (point.getY() >  ymax+EPS*um) || (point.getY() < ygap-EPS*um) )
	{
		streamlog_out(ERROR) << " position out of sensor! Y [mm] = "
			<< point.getY()/mm << " (Maximum y: " << ymax/mm 
			<< ", Minimum y:" << ygap << ")" << std::endl;
		isIn = false;
	}

	if( (point.getZ() > zmax+EPS*um) || (point.getZ() < -EPS*um) ) 
	{
		streamlog_out(ERROR) << " position out of sensor! Z [mm] = "
			<< point.getZ()/mm <<  " (Maximum z: " << zmax/mm << ")" << std::endl;
		isIn = false;
	}
	
	// Return if out or not
	return isIn;
}


// OTHER METHODS - LOCAL REF. SYSTEM

//
// Get info whether the given point is out of Si sensor (parameters: layerID,
// space point in local ref. system)
//
bool SiStripGeomFTD::isPointOutOfSensor(short int layerID, const CLHEP::Hep3Vector & point) const
{
	
	bool isOut = false;
	
	// Boundary set +- epsilon
	// Forward-type sensors
	double tanAlpha          = (getSensorWidth(layerID) - getSensorWidth2(layerID))/2./getSensorLength(layerID);
     	double actualSensorWidth = (getSensorWidth(layerID) - 2*tanAlpha*point.getZ());
      
        // Recalculate point into "local" local system (width depends on posZ)
    	CLHEP::Hep3Vector recalcPoint(point.getX(), point.getY() - 
    	  	      getSensorWidth(layerID)/2. + actualSensorWidth/2., point.getZ());
    	
    	if( (recalcPoint.getX() > (getSensorThick(layerID) +EPS*um)) || 
    	  	      (recalcPoint.getX() < (-EPS*um)) ||
    	  	      (recalcPoint.getY() > (actualSensorWidth+EPS*um)) 
    	  	      || (recalcPoint.getY() < (-EPS*um)) ||
    	  	      (recalcPoint.getZ() > (getSensorLength(layerID)+EPS*um)) ||
    	  	      (recalcPoint.getZ() < (-EPS*um)) )
    	{
    	        isOut = true;
    	}

    	// Return if out or not
    	return isOut;
}


//
// Get number of strips in sensors
//
int SiStripGeomFTD::getSensorNStrips(const int & diskID, const int & sensorID) const
{
	std::vector<int> const * sensorNStrips = 0;

	if(sensorID == 1 || sensorID == 2)
	{
		sensorNStrips = &_sensorNStripsInFront;
	}
	else if(sensorID == 3 || sensorID == 4)
	{
		sensorNStrips = &_sensorNStripsInRear;
	}
	else
	{
		streamlog_out(ERROR) << "SiStripGeomFTD::getSensorNStrips "
			<< " - Inconsistent ID of sensor: '" << sensorID << "'" << std::endl;
		exit(-1);
	}

	if(sensorNStrips->size()>(unsigned short int)diskID)
	{
		return (*sensorNStrips)[diskID];
	}
	else
	{
		return 0;
	}
}

//
// Get strip ID for RPhi  (perpendicular to the beam axis), 
// points are given in local ref. system.
//
//NOTES: Necesitem la posRPhi per designar quin  strip estic utilitzant,
//       la pos Z es necesaria per definir a quina zona del trapezoide 
//       estem, recorda que la Pitch= Pitch(Z)
/*int SiStripGeomFTD::getStripIDInRPhi(short int diskID, double posRPhi, double posZ ) const
{
	// Get pitch
	double sensPitch = getSensorPitchInRPhi(diskID,posZ);
	if(sensPitch == 0) 
	{
		streamlog_out(ERROR) << "SiStripGeomFTD::getStripIDInRPhi " 
			<< "- division by zero (sensPitch is zero)!!!"
			<< std::endl;
		exit(-1);
	}
	// Get number of strips
	int sensNStrips = getSensorNStripsInRPhi(diskID);
	
   	int stripID;
	// Calculate stripID
	if(posRPhi <= 0.)
	{
		stripID = 0;
	}
	else 
	{
		stripID = floor(posRPhi/sensPitch);
		if(stripID >= sensNStrips) 
		{
			stripID = sensNStrips - 1;
		}
	}
	// Error
	if(stripID >= sensNStrips) 
	{
		streamlog_out(ERROR) << "SiStripGeom::getStripIDInRPhi " 
			<< "- stripID in RPhi greater than number of strips!!!"
			<< std::endl;
		exit(-1);
	}
std::cout << " SiStripGeomFTD::getStripIDInRPhi" << std::endl;
std::cout << "----> " << posRPhi << std::endl;
std::cout << "----> " << sensPitch << std::endl;
std::cout << "----> " << sensNStrips << std::endl;
std::cout << "----> " << stripID << " (en zlocal=" << posZ << ")" << std::endl;
std::cout << "----> " << sensPitch*stripID << " (Y-position)" << std::endl;
std::cout << "END SiStripGeomFTD::getStripIDInRPhi" << std::endl;
	return stripID;
}*/

//
// Get the stereo angle given a sensor
// 
double SiStripGeomFTD::getStereoAngle(const int & diskID, const int & sensorID) const
{
	// Stereo angle (if proceed)
	double stAngle = 0.0;
	
	if(sensorID == 3 || sensorID == 4)
	{
		stAngle = 5e-3;
	}

	return stAngle;
}


//
// Get strip ID, point is given in local ref. system
// 
//
int SiStripGeomFTD::getStripID(const int & diskID, const int & sensorID,
		const double & posRPhi, const double & posZ) const
{
	//FIXME: BORRA
	double posRPhiRot=0.;
	double yorigen=0.0;
	// Get pitch
	double sensPitch = getSensorPitch(diskID,sensorID,posZ);
	if(sensPitch < 1e-40) 
	{
		streamlog_out(ERROR) << "SiStripGeomFTD::getStripID " 
			<< "- division by zero (sensPitch is zero)!!!"
			<< std::endl;
		exit(-1);
	}
	// Get number of strips
	int sensNStrips = getSensorNStrips(diskID,sensorID);
	
   	int stripID;
	// Calculate stripID
	if(posRPhi <= 0.)
	{
		stripID = 0;
	}
	else 
	{
		CLHEP::Hep3Vector pointRot = transformPointToRotatedLocal(diskID,sensorID,
				CLHEP::Hep3Vector(0.0,posRPhi,posZ) );
		posRPhiRot = pointRot.getY();
		// Put the strip 0 in the edge of the petal (recall
		// the y=0 of the local system begins in the long edge)
		const double tanPhi = tan(_ftdLayer->getPhiHalfDistance(diskID));
		yorigen = (getSensorLength(diskID)-pointRot.getZ())*tanPhi;

		stripID = floor((posRPhiRot-yorigen)/sensPitch);
		if(stripID >= sensNStrips) 
		{
			stripID = sensNStrips - 1;
		}
	}
	// Error
	if(stripID >= sensNStrips) 
	{
		streamlog_out(ERROR) << "SiStripGeom::getStripID " 
			<< "- stripID in RPhi greater than number of strips!!!"
			<< std::endl;
		exit(-1);
	}
std::cout << " SiStripGeomFTD::getStripID" << std::endl;
std::cout << "----> Y=" << posRPhi << " (y/Pitch):" << (posRPhi-yorigen)/sensPitch << std::endl;
std::cout << "----> YRot=" << posRPhiRot << " (y/Pitch):" << (posRPhiRot-yorigen)/sensPitch << std::endl;
std::cout << "----> Pitch=" << sensPitch << std::endl;
std::cout << "----> Number Strips:" << sensNStrips << std::endl;
std::cout << "----> StripID:" << stripID << " (en zlocal=" << posZ << ")" << std::endl;
std::cout << "----> NS*sID=" << yorigen+sensPitch*stripID << " (Y-position)" << std::endl;
std::cout << "END SiStripGeomFTD::getStripIDInRPhi" << std::endl;
	return stripID;
}

// METHODS TO IDENTIFY

//
// Get sensor pitch for RPhi strips (perpendicular to the beam axis).
// Currently it means the front face to the IP sensors of the petals.
// The posZ is needed because the trapezoid shape of the sensors therefore pitch=pitch(Z)
//
/*double SiStripGeomFTD::getSensorPitchInRPhi(short int diskID, double posZ) const
{
	if( (getLayerType(diskID) == strip) && 
			(_sensorPitchInFront.size() > (unsigned int)diskID) )
	{
		const double tanAlpha = tan(_ftdLayer->getPhiHalfDistance(diskID));
		if( (posZ*um>-EPS*um) && (posZ<(getSensorLength(diskID) + EPS*um)) ) //FIXME::::
		{
			return ((getSensorWidth(diskID) - 2.0*tanAlpha*posZ)
					/getSensorNStripsInRPhi(diskID));
		}
		else {
	std::cout << getSensorLength(diskID) << std::endl;
	std::cout << posZ << "   -EPS:" << -EPS*um << " Cumple?: " << (posZ*um>-EPS*um) << std::endl;
			streamlog_out(ERROR) 
				<< std::setiosflags(std::ios::fixed | std::ios::internal ) 
				<< std::setprecision(2)  
				<< "SiStripGeomFTD::getSensorPitchInRPhi - posZ: " 
				<< std::setw(4) << posZ << " out of range!!!" 
				<< std::resetiosflags(std::ios::showpos) 
				<< std::setprecision(0) 
				<< std::endl;
			exit(0);;
		}
	}
	else
	{
		return 0.;
	}
}*/

//
// Get sensor pitch
// The posZ is needed because the trapezoid shape of the sensors therefore pitch=pitch(Z)
// The posZ is in the local ref. system
//
double SiStripGeomFTD::getSensorPitch(const int & diskID, const int & sensorID, 
		const double & posZ) const
{
	if( (posZ < -EPS*um) && (posZ > (getSensorLength(diskID) + EPS*um)) )
	{
		streamlog_out(ERROR) 
			<< std::setiosflags(std::ios::fixed | std::ios::internal ) 
			<< std::setprecision(2)  
			<< "SiStripGeomFTD::getSensorPitch - posZ: " 
			<< std::setw(4) << posZ << " out of range!!!" 
			<< std::resetiosflags(std::ios::showpos) 
			<< std::setprecision(0) 
			<< std::endl;
		exit(0);
	}
	
	const double tanPhi = tan(_ftdLayer->getPhiHalfDistance(diskID));
//	CLHEP::Hep3Vector point(0.0,0.0,posZ);

/*	// Changing reference system: rotate the petal in its centroid point
	// to obtain the strips with the stereo angle
	double ycentre = getSensorWidth2(diskID)/2.0;
	if(getSensorWidth(diskID)/2.0 > ycentre)
	{
		ycentre = getSensorWidth(diskID)/2.0;
	}
	const double zcentre = getSensorLength(diskID)/2.0;

	// The centre of the rotation
	CLHEP::Hep3Vector rotcentre(0.0,ycentre,zcentre);
	// Translating Zpos to the new ref. system
	point -= rotcentre;

	// Rotating the posZ vector to the new coordinate system
	// Now we have the hit point from the system placed in the centre
	// of the sensor and rotated the stereo angle
	point.rotateX(-stAngle);
	// We need to return to the local ref. frame (defined along this code)
	// placed in a petal 'rotated'
	rotcentre.rotateX(-stAngle);
	// And extract the point in that system
	point += rotcentre;*/
	CLHEP::Hep3Vector point = transformPointToRotatedLocal( diskID, sensorID,
			CLHEP::Hep3Vector(0.,0.0,posZ) );

	// Now we can deal the petal in the usual way. 
	const double width_z = 2.0*(getSensorLength(diskID)-point.getZ())*tanPhi;
		
	// Note: to avoid incoherences
	double totalwidth = getSensorWidth2(diskID);
	if(getSensorWidth(diskID) > totalwidth)
	{
		totalwidth = getSensorWidth(diskID);
	}
std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
std::cout << " totalwidth= " << totalwidth << " width_z= " << width_z << " NStrips:" << getSensorNStrips(diskID,sensorID) << std::endl;
std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;

	return ((totalwidth-width_z)/getSensorNStrips(diskID,sensorID));
}

//
// Get the Y position (in the local ref. frame) for the given strip, posZ 
// is the point in the local reference system of the petal
//
double SiStripGeomFTD::getStripPos(const int & diskID, const int & sensorID, 
		const int & stripID, const double & posZ) const
{
	// Can't be stripID = 0
	if(stripID <= 0)
	{
		streamlog_out(ERROR) << "SoStripGeomFTD::getStripPos "
			<< "- incoherent stripID=" << stripID << "!!" 
			<< std::endl;
		exit(-1);
	}

	// Extract zPos position in the local rotated frame (if proceed)
	const double zPosPrime = transformPointToRotatedLocal( diskID,
			sensorID, CLHEP::Hep3Vector(0.0,0.0,posZ) ).getZ();

	// Get pitch
	double sensPitch = getSensorPitch(diskID,sensorID,zPosPrime);
	if(sensPitch < 1e-40) 
	{
		streamlog_out(ERROR) << "SiStripGeomFTD::getStripPos " 
			<< "- division by zero (sensPitch is zero)!!!"
			<< std::endl;
		exit(-1);
	}

	// Calculate position
	double posRPhi = sensPitch*(stripID + 0.5);
	// Error
	if ( (posRPhi<0.) || (posRPhi>getSensorWidth(diskID)) ) 
	{
		streamlog_out(ERROR) 
			<< "SiStripGeom::getStripPos - position out of sensor!!!"
			<< std::endl;
		exit(-1);
	}
	
	// Return R-Phi position of given strip in local ref. system
	return posRPhi;
}

//
// Transforming a given point to the local ref. frame of a petal which is rotated
// around its center an angle stAngle
//

CLHEP::Hep3Vector SiStripGeomFTD::transformPointToRotatedLocal(const int & diskID, 
		const int & sensorID, const CLHEP::Hep3Vector & point) const
{
	CLHEP::Hep3Vector pointPrime(point);

	const double stAngle = getStereoAngle(diskID,sensorID);

	// Changing reference system: rotate the petal in its centroid point
	// to obtain the strips with the stereo angle
	double ycentre = getSensorWidth2(diskID)/2.0;
	if(getSensorWidth(diskID)/2.0 > ycentre)
	{
		ycentre = getSensorWidth(diskID)/2.0;
	}
	const double zcentre = getSensorLength(diskID)/2.0;

	// The centre of the rotation
	CLHEP::Hep3Vector rotcentre(0.0,ycentre,zcentre);
	// Translating Zpos to the new ref. system
	pointPrime -= rotcentre;

	// Rotating the posZ vector to the new coordinate system
	// Now we have the hit point from the system placed in the centre
	// of the sensor and rotated the stereo angle
//	pointPrime.rotateX(-stAngle);
	// We need to return to the local ref. frame (defined along this code)
	// placed in a petal 'rotated'
	rotcentre.rotateX(-stAngle);
	// And extract the point in that system
	pointPrime += rotcentre;

	pointPrime.rotateX(-stAngle);
	/* // Consistency Check
	CLHEP::Hep3Vector Orig(0.,0.,0.);
	Orig = CLHEP::Hep3Vector(0.0,ycentre,zcentre)-rotcentre;
	std::cout << " pointPrime+O =? P :" << pointPrime + Orig << " =? " << point << std::endl;*/
	
	return pointPrime;
}


// PRINT METHODS

//
// Method printing general Gear parameters
//

void SiStripGeomFTD::printGearParams() const
{
   streamlog_out(MESSAGE3) << std::endl
                           << " "
	                   		<< DUNDERL
                           << DBLUE
                           << "Gear parameters:"
                           << ENDCOLOR
                           << " "
                           << std::endl  << std::endl;

   // Gear type: BField
   streamlog_out(MESSAGE3) << "  B field [T]:       "
                           << _magField/T
                           << std::endl << std::endl;

   // Gear type: VXD
   if (_gearType == "VXD") {

      // Print general info
      for (int i=0; i<_numberOfLayers; ++i){

         streamlog_out(MESSAGE2) << std::endl;
         streamlog_out(MESSAGE2) << "  Layer: "           << _layerRealID[i]      << std::endl;
         if (_layerType[i]==pixel) {
         streamlog_out(MESSAGE2) << "  LayerType:       " << "pixel"              << std::endl;
         }
         else {
         streamlog_out(MESSAGE2) << "  LayerType:       " << "strip"              << std::endl;
         }
         streamlog_out(MESSAGE2) << "  NumberOfLadders: " << _numberOfLadders[i]  << std::endl;
         streamlog_out(MESSAGE2) << "  Radius[mm]:      " << _layerRadius[i]/mm   << std::endl;
         streamlog_out(MESSAGE2) << "  Width[mm]:       " << _ladderWidth[i]/mm   << std::endl;
         streamlog_out(MESSAGE2) << "  Length[mm]:      " << _ladderLength[i]/mm  << std::endl;
         streamlog_out(MESSAGE2) << "  Petal SemiAngle: " << _layerHalfPhi[i]     << std::endl;
         streamlog_out(MESSAGE2) << "  Phi0:            " << _layerPhi0[i]        << std::endl;
         streamlog_out(MESSAGE2) << "  Theta:           " << _layerTheta[i]       << std::endl;
         streamlog_out(MESSAGE2) << "  OffsetY[mm]:     " << _ladderOffsetY[i]/mm << std::endl;
 //        streamlog_out(MESSAGE2) << "  OffsetZ[mm]:     " << _ladderOffsetZ[i]/mm << std::endl; OUT//
      }

   }
   // Gear type: unknown - error
   else {
   	streamlog_out(ERROR) << "Unknown gear type!"
   	                     << std::endl;

      exit(0);
   }
}
//
// Method printing sensor Gear parameters (parameters: sensorID)
//
void SiStripGeomFTD::printSensorParams(short int layerID) const
{
	// Print sensor parameters
	streamlog_out(MESSAGE2) << "  * Parameters: " << _sensorThick[layerID]/um       << "um thick, "
		<< "with "            << _sensorPitchInRear[layerID]/um    << "um pitch "
		<< "and "             << _sensorNStripsInRear[layerID]     << " strips in Z"
		<< ", resp. "         << _sensorPitchInFront[layerID]/um << "um pitch "
		<< "and "             << _sensorNStripsInFront[layerID]  << " strips in R-Phi."
		<< std::endl;
}

} // Namespace;

