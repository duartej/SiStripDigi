
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

// FIXME: Valores en cm o en mm????
//        Parece que en la base de datos lo he
//        guardado en mm y aqui utiliza cm como
//        unidad...


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
	_sensorPitchInRPhi  = std::vector<double>(2*_numberOfLayers,50.0 * um); //FIXME: HARDCODED
	_sensorNStripsInRPhi.reserve(2*_numberOfLayers);
	for(int i = 0; i < _numberOfLayers; i++)
	{
		const double xmaxsensorup = _ftdLayer->getSensitiveLengthMax(i)*mm;
		_sensorNStripsInRPhi.push_back( (int)(xmaxsensorup/_sensorPitchInRPhi[i])+1 );
	}
	_sensorNStripsInRPhi.insert(_sensorNStripsInRPhi.end(),_sensorNStripsInRPhi.begin(),
			_sensorNStripsInRPhi.end());

	_sensorPitchInZ     = std::vector<double>(2*_numberOfLayers,50.0*um); // FIXME
	_sensorNStripsInZ.reserve(2*_numberOfLayers);
	for(int i = 0; i < _numberOfLayers; i++)
	{
		const double widthsensor = _ftdLayer->getSensitiveWidth(i)*mm;
		_sensorNStripsInZ[i] = (int)(widthsensor/_sensorPitchInZ[i])+1;
	}
	_sensorNStripsInZ.insert(_sensorNStripsInZ.end(),_sensorNStripsInZ.begin(),
			_sensorNStripsInZ.end());
}

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

//
// Method transforming given point from global ref. system to local ref. system
// FIXME: Description of the local reference system---
// (parameters: diskID, petalID, sensorID and space point in global ref. system)
//
CLHEP::Hep3Vector SiStripGeomFTD::transformPointToLocal(short int diskID, short int petalID, short int sensorID, const CLHEP::Hep3Vector & globalPoint)
{
//FIXME: PROV
	// Initialize local point
	CLHEP::Hep3Vector localPoint(globalPoint);

std::cout << " transformPointToLocal: processing HIT" << std::endl;
std::cout << " --- No processed hit-point: " << localPoint<< std::endl;
	// Phi Angle of the Petal (on the centroid of the petal with respect the x-axis)
	const double phi0 = _ftdLayer->getPhiPetalCd(diskID,petalID);
std::cout << "Phi petal: " << phi0*180.0/M_PI << std::endl;
	
	//
        // The local ref. frame (LRF') is defined positive and with
        // the unitary vectors as
        //    x'=direction of the electric field of the IP faced sensors 
        //       (for all faced sensors, the global Z-axis)
        //    y'=paralel to the Z-Strips (paralel to the smallest dimension
        //       of the petal)
        //    z'=paralel to the Phi-Strips (paralel to the largest dimension
        //       of the petal)
	// Local ref. system it is defined in a way as the X,Y,Z 
	// of a hit is always positive in this frame. So, the 
	// (0,0,0) in this local frame correspond to the point
	// remaining in the the trapezoid's smallest side and
        // in the back face w.r.t. the IP.
        //
	
	// 
	// RECALL: ftd gear have the convention mm (distance)
	// 

	// Extract the Z of the sensor (in the back face)
	const double zsensor = _ftdLayer->getSensitiveZposition(diskID,petalID,sensorID)*mm;
	int zsign = (int)(zsensor/fabs(zsensor));
	const double sensorthickness = _ftdLayer->getSensitiveThickness(diskID)*mm;
	const double zsensorback = zsensor+zsign*sensorthickness/2.0;
	// And the X and Y CentroiD position: over the smallest side of the trapezoid
	const double xsensorCd = (_ftdLayer->getSensitiveRinner(diskID)*mm)*cos(phi0);
	const double ysensorCd = (_ftdLayer->getSensitiveRinner(diskID)*mm)*sin(phi0);
std::cout << "\n --- Sensor centroid point: (" << xsensorCd << "," << ysensorCd <<
	 "," << zsensorback << ")" ;
std::cout << "  diskID:" << diskID << "(Real ID:" 
	<< _layerRealID.at(diskID) << ") petalID:" << 
	petalID << " sensorID:" << sensorID <<	std::endl;

	// The (0,0,0) position of LFR
	CLHEP::Hep3Vector localOrigin(xsensorCd, ysensorCd, zsensorback);

	// Translating the globalPoint to the local Origin
	localPoint -= localOrigin;

std::cout << " --- Antes: " << localPoint;
	
	// Perform rotation to get the system local
	localPoint.rotateZ(-zsign*phi0);
	localPoint.rotateY(-zsign*M_PI/2.0);
std::cout << " --- Despues: " << localPoint;

	// Avoiding X,Y and Z negative values--> displacing from the
        // centroid to the edge
        const double longtrapezoidedge = _ftdLayer->getSensitiveLengthMax(diskID)*mm;
	localPoint += CLHEP::Hep3Vector(0.0,longtrapezoidedge/2.0,0.0); 
std::cout << " --- TOTAL: " << localPoint;
	
	// Return space point in local ref. system
	return localPoint;
}

//
// Method transforming given vector from global ref. system to local ref. system
// (parameters: layerID, ladderID, sensorID and vector in global ref. system)
//
CLHEP::Hep3Vector SiStripGeomFTD::transformVecToLocal(short int diskID, short int petalID, const CLHEP::Hep3Vector & globalVec)
{
//FIXME: PROV
	// Initialize local vector
	CLHEP::Hep3Vector localVec(globalVec);

std::cout << " VECTOR TRANS: " << localVec << " ---> ";
	
	const double phi0 = _ftdLayer->getPhiPetalCd(diskID,petalID);
	// Perform rotation to get the system local
	localVec.rotateZ(-phi0);
	localVec.rotateY(-M_PI/2.0);
	
	
	// Return vector in local ref. system
std::cout << localVec << std::endl;
	return localVec;
}

//
// Method transforming given matrix 3x3 from global ref. system to local ref. system
// (parameters: layerID, ladderID, sensorID and matrix in global ref. system)
//
CLHEP::HepMatrix SiStripGeomFTD::transformMatxToLocal(short int layerID, short int ladderID, const CLHEP::HepMatrix & globalMatrix)
{
   // Initialize local matrix 3x3 to zero values
	CLHEP::HepMatrix localMatrix(3,3,0);

   // Initialize rotation matrices: R, R^T (transposition)
	CLHEP::HepMatrix rotMatrix(3,3,0);
	CLHEP::HepMatrix rotMatrixT(3,3,0);

	// Gear type: VXD
	   if (_gearType == "VXD") {

	      // Calculate rotation angles
	      double theta = getLadderTheta(layerID);
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
CLHEP::Hep3Vector SiStripGeomFTD::transformPointToGlobal(short int layerID, short int ladderID, short int sensorID, const CLHEP::Hep3Vector & localPoint)
{
   // Initialize global point
   CLHEP::Hep3Vector globalPoint(localPoint);

   // Gear type: VXD
   if (_gearType == "VXD") {

      // Calculate rotation angles
      double theta = getLadderTheta(layerID);
      double phi   = getLadderPhi(layerID, ladderID);

      // Find (0,0,0) position of local coordinate system
      CLHEP::Hep3Vector localOrigin(getLayerRadius(layerID), getLadderOffsetY(layerID), getLadderOffsetZ(layerID));
      localOrigin.rotateZ(+phi);

      // Perform translation - to the center of a ladder
      globalPoint -= CLHEP::Hep3Vector(+getSensorThick(layerID)/2., +getSensorWidth(layerID)/2.,
                                +0.5*getLadderLength(layerID) - (2*sensorID + 1)*getSensorRimWidthInZ(layerID) - sensorID*getSensorGapInBetween(layerID)
                                - sensorID*getSensorLength(layerID));

      // Perform rotation - with respect to the local centre of a ladder
      globalPoint.rotateY(-theta);
      globalPoint.rotateZ(+phi);

      // Perform translation - to the global system
      globalPoint += localOrigin;

   }

   // Gear type: unknown - error
   else {
      streamlog_out(ERROR) << "Unknown gear type!"
                           << std::endl;

      exit(0);
   }

   // Return space point in global ref. system
	return globalPoint;
}

//
// Method transforming given vector from local ref. system to global ref. system
// (parameters: layerID, ladderID, sensorID and vector in local ref. system)
//
CLHEP::Hep3Vector SiStripGeomFTD::transformVecToGlobal(short int layerID, short int ladderID, const CLHEP::Hep3Vector & localVec)
{
   // Initialize global vector
	CLHEP::Hep3Vector globalVec(localVec);

   // Gear type: VXD
   if (_gearType == "VXD") {

      // Calculate rotation angles
      double theta = getLadderTheta(layerID);
      double phi   = getLadderPhi(layerID, ladderID);

      // Perform rotation - to the center of a ladder
      globalVec.rotateY(-theta);
      globalVec.rotateZ(+phi);

   }

   // Gear type: unknown - error
   else {
      streamlog_out(ERROR) << "Unknown gear type!"
                           << std::endl;

      exit(0);
   }

   // Return vector in global ref. system
   return globalVec;
}

//
// Method transforming given matrix 3x3 from local ref. system to global ref. system
// (parameters: layerID, ladderID, sensorID and matrix in local ref. system)
//
CLHEP::HepMatrix SiStripGeomFTD::transformMatxToGlobal(short int layerID, short int ladderID, const CLHEP::HepMatrix & localMatrix)
{
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
bool SiStripGeomFTD::isPointInsideSensor (short int layerID, short int ladderID, short int sensorID, const CLHEP::Hep3Vector & point) const
{
   bool isIn = false;

   // Gear type: VXD
   if (_gearType == "VXD") {
      //
      // Transform point to local ref. system

      // Initialize local point
      CLHEP::Hep3Vector localPoint(point);

      // Calculate rotation angles
      double theta = getLadderTheta(layerID);
      double phi   = getLadderPhi(layerID, ladderID);

      // Find (0,0,0) position of local coordinate system
      CLHEP::Hep3Vector localOrigin(getLayerRadius(layerID), getLadderOffsetY(layerID), getLadderOffsetZ(layerID));
      localOrigin.rotateZ(+phi);

      // Perform translation - to the center of a ladder
      localPoint -= localOrigin;

      // Perform rotation - to the center of a ladder
      localPoint.rotateZ(-phi);
      localPoint.rotateY(+theta);

      // Perform translation such as X, Y, Z are positive
      localPoint += CLHEP::Hep3Vector(+getSensorThick(layerID)/2., +getSensorWidth(layerID)/2.,
                               +0.5*getLadderLength(layerID) - (2*sensorID + 1)*getSensorRimWidthInZ(layerID) - sensorID*getSensorGapInBetween(layerID)
                               - sensorID*getSensorLength(layerID));



      // Boundary set +- epsilon
      // Barrel-type sensors
      if ( (getLayerType(layerID) == stripB) && ( (localPoint.getX() < (getSensorThick(layerID) +EPS*um)) && (localPoint.getX() > (-EPS*um)) &&
                                                  (localPoint.getY() < (getSensorWidth(layerID) +EPS*um)) && (localPoint.getY() > (-EPS*um)) &&
                                                  (localPoint.getZ() < (getSensorLength(layerID)+EPS*um)) && (localPoint.getZ() > (-EPS*um)) ) ) isIn = true;
      // Forward-type sensors
      if (getLayerType(layerID) == stripF) {

         double tanAlpha          = (getSensorWidth(layerID) - getSensorWidth2(layerID))/2./getSensorLength(layerID);
         double actualSensorWidth = (getSensorWidth(layerID) - 2*tanAlpha*localPoint.getZ());

         // Recalculate point into "local" local system (width depends on posZ)
         CLHEP::Hep3Vector recalcPoint(localPoint.getX(), localPoint.getY() - getSensorWidth(layerID)/2. + actualSensorWidth/2., localPoint.getZ());

         if ( (recalcPoint.getX() < (getSensorThick(layerID) +EPS*um)) && (recalcPoint.getX() > (-EPS*um)) &&
              (recalcPoint.getY() < (actualSensorWidth       +EPS*um)) && (recalcPoint.getY() > (-EPS*um)) &&
              (recalcPoint.getZ() < (getSensorLength(layerID)+EPS*um)) && (recalcPoint.getZ() > (-EPS*um)) ) isIn = true;
      }

   }
   // Gear type: unknown - error
   else {
      streamlog_out(ERROR) << "Unknown gear type!"
                           << std::endl;

      exit(0);
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
// Get strip ID for strips perpendicular to the RPhi local ref. frame plane, 
// points are given in local ref. system.
//
//NOTES: Necesitem la posRPhi per designar quin  strip estic utilitzant,
//       la pos Z es necesaria per definir a quina zona del trapezoide 
//       estem, recorda que la Nstrips= Nstrips(Z)
int SiStripGeomFTD::getStripIDInRPhi(short int diskID, double posRPhi, double posZ ) const
{
	// Get pitch
	double sensPitch = getSensorPitchInRPhi(diskID);
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
/*std::cout << " SiStripGeomFTD::getStripIDInRPhi" << std::endl;
std::cout << "----> " << sensPitch << std::endl;
std::cout << "----> " << sensNStrips << std::endl;
std::cout << "----> " << stripID << " (en zlocal=" << posZ << ")" << std::endl;
std::cout << "END SiStripGeomFTD::getStripIDInRPhi" << std::endl;*/
	return stripID;
}

// METHODS TO IDENTIFY

//
// Get sensor pitch for strips perpendicular to the RPhi plane in the local ref. system.
// Currently it means the faced to IP sensors of the petals.
// The posZ is needed because the trapezoid shape of the sensors therefore pitch=pitch(Z)
//
double SiStripGeomFTD::getSensorPitchInRPhi(short int diskID, double posZ) const
{
	if( (getLayerType(diskID) == strip) && 
			(_sensorPitchInRPhi.size() > (unsigned int)diskID) )
	{
		const double tanAlpha = tan(_ftdLayer->getPhiHalfDistance(diskID));
		if( (posZ>-EPS*um) && (posZ<(getSensorLength(diskID) + EPS*um)) ) 
		{
			return ((getSensorWidth(diskID) - 2.0*tanAlpha*posZ)
					/getSensorNStripsInRPhi(diskID));
		}
		else {
			streamlog_out(ERROR) 
				<< std::setiosflags(std::ios::fixed | std::ios::internal ) 
				<< std::setprecision(2)  
				<< "SiStripGeom::getSensorPitchInRPhi - posZ: " 
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
		<< "with "            << _sensorPitchInZ[layerID]/um    << "um pitch "
		<< "and "             << _sensorNStripsInZ[layerID]     << " strips in Z"
		<< ", resp. "         << _sensorPitchInRPhi[layerID]/um << "um pitch "
		<< "and "             << _sensorNStripsInRPhi[layerID]  << " strips in R-Phi."
		<< std::endl;
}

} // Namespace;

