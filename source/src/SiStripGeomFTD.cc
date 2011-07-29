
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
#include <gear/GearParameters.h>
//#include <gear/FTDParameters.h>
//#include <gear/FTDLayerLayout.h>
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
	SiStripGeom(detector)
{
}

SiStripGeomFTD::~SiStripGeomFTD()
{
}


//FIXME
//FIXME: PROVISIONAL HASTA QUE IMPLEMENTE LOS FTD EN GEAR
//FIXME

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
	catch (gear::UnknownParameterException& e) 
	{
		std::cout << "No magnetic field found in gear file!" << std::endl;
	}

	//------Get the geometry from the gear file-----// FIXME: PROVISIONAL HASTA que se Implemente FTD en gear
	//FIXME!!! PROVISIONAL hasta que se implemente FTD en gear
	const gear::GearParameters & ftdgearParams = marlin::Global::GEAR->getGearParameters(this->_gearType);

	_numberOfLayers  = ftdgearParams.getIntVals("DiskNumbers").size();
	_layerZ          = ftdgearParams.getDoubleVals("FTDZCoordinate");
	_layerRadius     = ftdgearParams.getDoubleVals("FTDInnerRadius");
	_layerOuterRadius= ftdgearParams.getDoubleVals("FTDOuterRadius");
	_layerRealID     = ftdgearParams.getIntVals("DiskNumbers"); //Posiblemente inecesario, a partir de FTDZCoordinatea
	// From mm to cm
	convertunits(_layerZ,(double)mm);
	convertunits(_layerRadius,(double)mm);
	convertunits(_layerOuterRadius,(double)mm);
	//
	_numberOfLadders= ftdgearParams.getIntVals("NumberOfLadders");

	_sensorThick        = ftdgearParams.getDoubleVals("FTDDiskSiThickness");
	_sensorPitchInZ     = std::vector<double>(_numberOfLayers,-1); // NO HAY EN Z
	_sensorNStripsInZ   = std::vector<int>(_numberOfLayers,-1); // NO HAY EN Z
	_sensorPitchInRPhi  = std::vector<double>(_numberOfLayers,50.0 * um); //FIXME: HARDCODED
	_sensorWidth        = ftdgearParams.getDoubleVals("UpSensordxMax");  // x-direction
	_sensorLength       = ftdgearParams.getDoubleVals("FTDPetalDy");     // y-direction

	_ladderOffsetY      = ftdgearParams.getDoubleVals("SensorOffsetY");  // Ya es la mitad de 5
	_ladderOffsetZ      = ftdgearParams.getDoubleVals("SensorThickness");
	convertunits(_ladderOffsetY,(double)mm);
	convertunits(_ladderOffsetZ,(double)mm);
	convertunits(_sensorWidth,(double)mm);  
	convertunits(_sensorThick,(double)mm);
	convertunits(_sensorLength,(double)mm);
	
	_ladderLength       = _sensorLength;
	// FIXME PROVISONAL: 
	for(unsigned int i=0; i < _ladderOffsetZ.size(); i++)
	{
		_ladderOffsetZ[i] = _ladderOffsetZ[i]/2.0;
	}

	std::vector<double> xminsensordown = ftdgearParams.getDoubleVals("DownSensordxMin");
	convertunits(xminsensordown,(double)mm);
	for(unsigned int i = 0; i < xminsensordown.size(); i++)
	{
		_sensorNStripsInRPhi.push_back( (int)(xminsensordown[i]/_sensorPitchInRPhi[i])+1 );
	}

	_layerPhi0   = ftdgearParams.getDoubleVals("FTDPetalPhi0");
	const double degree = M_PI/180.0;
	convertunits(_layerPhi0,degree);
}

// FIXME: Provisional mientras no esten codificada la FTD, pues son funciones lentas
std::map<std::string,short int> SiStripGeomFTD::cellIDDecProv(EVENT::SimTrackerHit * & simHit)
{
	const double hitpos[3] = { *(simHit->getPosition())*mm, *(simHit->getPosition()+1)*mm,
		*(simHit->getPosition()+2)*mm };
	std::map<std::string,short int> cellid;
	//Disk (layer)
	//FIXME: Hardcoded some extreme thickness (position of the sensors)
	double thickness = 5.0*mm;
	for(unsigned int i = 0; i < _layerZ.size(); i++)
	{
		short int sign = 1;
		if( hitpos[2] < 0.0 )
		{
			sign = -1;
		}
		double diff = sign*hitpos[2] - _layerZ[i];
		
		if( diff > thickness )
		{
			continue;
		}
		
		cellid["layer"] = sign*(i); //+1); FIXME!!
		// FIXME: Estimation!!
		// Ladder (Petal): increasing with theta
		const double thetahit   = atan(hitpos[1]/hitpos[0]);
		const double thetaPetals=  2.0*M_PI/_numberOfLadders[i];
		double idx = (thetahit-thetaPetals/2.0)/thetaPetals;
		short int retid;
		if( idx < 0.0 )
		{
			retid = int(idx) + 18;  // FIXME HARDCODED!!
			if( idx > -1 )
			{
				retid = 1;
			}
		}
		else
		{
			retid = int(idx)+1;
		}
		cellid["ladder"] = retid;
		// FIXME: Estimation!!
		// Sensor UP-DOWN and FRONT-REAR
		// cara al IP y de arriba a abajo: 1,2,3,4

		// La mitad
		const double rhalf = _layerRadius[i]+(_layerOuterRadius[i]-_layerRadius[i])/2.0;
		const double rhit  = sqrt(hitpos[0]*hitpos[0]+hitpos[1]*hitpos[1]);
		// UP or DOWN
		int idxSensor = 1; // Up- per default: multiple de 1
		if( rhit < rhalf )
		{
			idxSensor = 2;  //DOwN : multiple de dos
		}
		// FRONT (mirando al IP) or REAR, por defecto FRONT
		if( fabs(hitpos[2]) > _layerZ[i] )
		{
			idxSensor += 2;
		}
		cellid["sensor"] = idxSensor;


		break;
	}
	
	return cellid;
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
// (parameters: layerID, ladderID, sensorID and space point in global ref. system)
//
CLHEP::Hep3Vector SiStripGeomFTD::transformPointToLocal(short int layerID, short int ladderID, short int sensorID, const CLHEP::Hep3Vector & globalPoint)
{
	// Initialize local point
	CLHEP::Hep3Vector localPoint(globalPoint);
std::cout << "ANtes de nada: " << localPoint << std::endl;
	
      	// Calculate rotation angles
	double theta = 0.0; //getLadderTheta(layerID);
	double phi   = getLadderPhi(layerID, ladderID-1);

	// Find (0,0,0) position of local coordinate system
	// Para mantener que siempre tengamos x,y,z positivos en el sis.ref.local del sensor
        // utilizo el dxMax
	//CLHEP::Hep3Vector localOrigin(getLayerRadius(layerID), getLadderOffsetY(layerID), getLadderOffsetZ(layerID));
	double radius = getLayerRadius(layerID);
	double zpos   = _layerZ[layerID];
	if( sensorID < 3 )
	{
		zpos-= getLadderOffsetZ(layerID);
	}
	else
	{
		zpos+= getLadderOffsetZ(layerID);
	}
	CLHEP::Hep3Vector localOrigin(radius+getLadderOffsetX(layerID), 
			radius+getLadderOffsetY(layerID), zpos);
std::cout << " Origen del petalo: " << localOrigin << std::endl;
	localOrigin.rotateZ(+phi);
std::cout << " Origen del petalo (rotado +phi around z): " << localOrigin << std::endl;

	// Perform translation - to the center of a ladder
	localPoint -= localOrigin;
std::cout << " Translacion (-Origen del petalo): " << localPoint << std::endl;
	
	// Perform rotation - to the center of a ladder
	localPoint.rotateZ(-phi);
	localPoint.rotateY(+theta);

std::cout << " Rotacion (-phi) :" << localPoint << std::endl;
	
	// Perform translation such as X, Y, Z are positive
	//localPoint += CLHEP::Hep3Vector(+getSensorThick(layerID)/2., +getSensorWidth(layerID)/2.,
	//		+0.5*getLadderLength(layerID) - (2*sensorID + 1)*getSensorRimWidthInZ(layerID) 
	//		- sensorID*getSensorGapInBetween(layerID) - sensorID*getSensorLength(layerID));
	localPoint += CLHEP::Hep3Vector(+getLadderOffsetX(layerID)/2., +getLadderLength(layerID)/2.,
			+2.0*getSensorThick(layerID));
std::cout << " Final, donde xyz positvas: " << localPoint<< std::endl;
std::cout << " UN DESASTRE!!! " << localPoint<< std::endl;

      	// Check if local point within sensor boundaries +- epsilon
      if (isPointOutOfSensor(layerID, localPoint)) {

         streamlog_out(ERROR) << std::setprecision(3) << "SiStripGeomFTD::transformPointToLocal - point: "
                              << localPoint           << " is out of sensor!!!"
                              << std::setprecision(0) << std::endl;
         exit(0);
      }
   
   // Return space point in local ref. system
   return localPoint;
}

//
// Method transforming given vector from global ref. system to local ref. system
// (parameters: layerID, ladderID, sensorID and vector in global ref. system)
//
CLHEP::Hep3Vector SiStripGeomFTD::transformVecToLocal(short int layerID, short int ladderID, const CLHEP::Hep3Vector & globalVec)
{
   // Initialize local vector
	CLHEP::Hep3Vector localVec(globalVec);

   // Gear type: VXD
   if (_gearType == "VXD") {

      // Calculate rotation angles
      double theta = getLadderTheta(layerID);
      double phi   = getLadderPhi(layerID, ladderID);

      // Perform rotation - to the center of a ladder
      localVec.rotateZ(-phi);
      localVec.rotateY(+theta);

   }

   // Gear type: unknown - error
   else {
      streamlog_out(ERROR) << "Unknown gear type!"
                           << std::endl;

      exit(0);
   }

   // Return vector in local ref. system
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

