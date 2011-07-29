#include "SiStripGeom.h"

// Include standard header files
#include "Colours.h"
#include "PhysicalConstants.h"

#include <cstdlib>
#include <iomanip>
#include <algorithm>

// Include CLHEP header files
#include <CLHEP/Vector/EulerAngles.h>
#include <CLHEP/Vector/Rotation.h>

// Include Gear header files
//#include <gear/BField.h>
//#include <gear/VXDParameters.h>
//#include <gear/VXDLayerLayout.h>
#include <gearimpl/Vector3D.h>

// Include Marlin
#include <marlin/Global.h>
#include <streamlog/streamlog.h>

// Namespaces
//using namespace CLHEP;
//using namespace marlin;

namespace sistrip 
{
	SiStripGeom::SiStripGeom(const std::string & subdetector)
		:_gearType(subdetector)
	{

	}
	
	SiStripGeom::~SiStripGeom()
	{
	//	std::cout << "Deleting SiStripGeom" << std::endl;
	}
	
/* Virtual pure	//
	// Method initializing this class - reads Gear parameters from XML file
	//
	void SiStripGeom::initGearParams()
	{
		//------Get the geometry from the gear file-----//
		const gear::VXDParameters& gearVXD = Global::GEAR->getVXDParameters() ;
		const gear::VXDLayerLayout& layerVXD = gearVXD.getVXDLayerLayout(); 
	
	//number of layers
	_numberOfLayers  = layerVXD.getNLayers();
	
	//the number of ladders within layers
	_laddersInLayer.resize(_numberOfLayers);

	//Azimuthal offset of the whole structure of each layer.
	_layerHalfPhi.resize(_numberOfLayers);
	
	//layer thickness and half-thickness
	_layerHalfThickness.resize(_numberOfLayers);
	_layerThickness.resize(_numberOfLayers);
	
	//Distance from the middle of the Si sensor to IP in each layer.
	_layerRadius.resize(_numberOfLayers);
	
	//The length of the Si sensor in each layer
	_layerLadderLength.resize(_numberOfLayers);
	
	//The width and half-width of the Si sensor in each layer
	_layerLadderWidth.resize(_numberOfLayers);
	_layerLadderHalfWidth.resize(_numberOfLayers);

	//The offset of the sensitive area in ladder in each layer
	_layerActiveSiOffset.resize(_numberOfLayers);
	
	// The gaps in z between two subladders within each layer
	const gear::GearParameters& gearVXDInfra = Global::GEAR->getGearParameters("VXDInfra") ;
	const std::vector<double> laddergaps = gearVXDInfra.getDoubleVals("LadderGaps");
	std::cout<<"laddergaps size "<<laddergaps.size()<<std::endl;
	_layerLadderGap.resize(laddergaps.size());

	//ladder offset in phi
	_layerPhiOffset.resize(_numberOfLayers);
	
	for(int layer = 0; layer < _numberOfLayers; layer++)
	{
		_laddersInLayer[layer] = layerVXD.getNLadders(layer);  
		_layerHalfPhi[layer] = layerVXD.getPhi0(layer); 
		
		//should half phi be replaced by this from Alexei's original code?
		//_layerHalfPhi[layer] = PI/((double)_laddersInLayer[layer]); OJOOO
		
		_layerThickness[layer] =layerVXD.getSensitiveThickness(layer);
		_layerHalfThickness[layer] = 0.5*_layerThickness[layer];      
		_layerRadius[layer] =layerVXD.getSensitiveDistance(layer) + 0.5 * _layerThickness[layer];
		_layerLadderLength[layer] = 2*layerVXD.getSensitiveLength(layer);
		_layerLadderWidth[layer] = layerVXD.getSensitiveWidth(layer);
		_layerLadderHalfWidth[layer] = _layerLadderWidth[layer]/2.;
		_layerActiveSiOffset[layer] = - (layerVXD.getSensitiveOffset(layer));
		_layerLadderGap[layer] = laddergaps[layer];
		_layerPhiOffset[layer] = layerVXD.getPhi0(layer);
	}
		
	std::cout<<" _numberOfLayers "<<_numberOfLayers<<std::endl;
	std::cout<<" _pixelSizeX "<<_pixelSizeX<<std::endl;
	std::cout<<" _pixelSizeY "<<_pixelSizeY<<std::endl;
	std::cout<<" _electronsPerKeV "<<_electronsPerKeV<<std::endl;
	std::cout<<" _segmentDepth "<<_segmentDepth<<std::endl;
	std::cout<<" _currentTotalCharge "<<_currentTotalCharge<<std::endl;
	for (int i=0; i<_numberOfLayers; ++i) 
	{
		std::cout<<"layer "<<i<<std::endl;
		std::cout<<" _laddersInLayer "<<_laddersInLayer[i]<<std::endl;
		std::cout<<" _layerRadius "<<_layerRadius[i]<<std::endl;
		std::cout<<" _layerLadderLength "<<_layerLadderLength[i]<<std::endl;
		std::cout<<" _layerLadderHalfWidth "<<_layerLadderHalfWidth[i]<<std::endl;
		std::cout<<" _layerPhiOffset "<<_layerPhiOffset[i]<<std::endl;
		std::cout<<" _layerActiveSiOffset "<< _layerActiveSiOffset[i]<<std::endl;
		std::cout<<" _layerHalfPhi "<<_layerHalfPhi[i]<<std::endl;
		std::cout<<" _layerLadderGap "<<_layerLadderGap[i]<<std::endl;
		std::cout<<" _bkgdHitsInLayer "<<_bkgdHitsInLayer[i]<<std::endl;
		std::cout<<" _layerLadderWidth "<<_layerLadderWidth[i]<<std::endl;
		std::cout<<" _layerThickness "<<_layerThickness[i]<<std::endl;
		
		std::cout<<" _layerHalfThickness "<<_layerHalfThickness[i]<<std::endl;
	}
	
	//SCALING = 25000.;
}*/ 


// GEOMETRY PROPERTIES

//
// Encode cell ID
//
int SiStripGeom::encodeCellID(short int layerID, short int ladderID, short int sensorID) const
{
	return layerID*LAYERCOD + ladderID*LADDERCOD + sensorID*SENSORCOD;
}

//
// Decode cell ID
//
void SiStripGeom::decodeCellID(short int & layerID, short int & ladderID, short int & sensorID, int cellID) const
{
	layerID  =  cellID / LAYERCOD;
	ladderID = (cellID - layerID*LAYERCOD) / LADDERCOD;
	sensorID = (cellID - layerID*LAYERCOD - ladderID*LADDERCOD) / SENSORCOD;
}

//
// Encode strip ID
//
int SiStripGeom::encodeStripID(StripType type, int stripID) const
{
	// TODO: Posible mejora del algoritmo
	int encodedStripID = 0;
	
	// Strip in R-Phi
	if (type == RPhi) 
	{
		encodedStripID = (stripID + STRIPOFF)*STRIPCODRPHI;
	}
	
	// Strip in Z
	if (type == Z) 
	{
		encodedStripID = (stripID + STRIPOFF)*STRIPCODZ;
	}

	return encodedStripID;
}

//
// Decode strip ID
//
void SiStripGeom::decodeStripID(StripType & type, int & stripID, int encodedStripID) const 
{
	// TODO: Posible mejora del algortimo...

	// Strip in RPhi
	if (encodedStripID/STRIPCODRPHI > 0.) 
	{
		stripID = encodedStripID/STRIPCODRPHI - STRIPOFF;
		type    = RPhi;
		
		return;
	}
	// Strip in Z
	if (encodedStripID/STRIPCODZ > 0.) 
	{
		stripID = encodedStripID/STRIPCODZ - STRIPOFF;
		type    = Z;
		
		return;
	}

	// Error
	streamlog_out(ERROR) << "SiStripGeom::decodeStripID: "
		<< encodedStripID
		<< " - problem to identify if strips in Z or R-Phi!!!"
		<< std::endl;
	exit(0);
}

//
// Get layer real ID
//
int SiStripGeom::getLayerRealID(short int layerID) const
{
	//TODO: Posible mejora del algoritmo: try-catch
	//      capturando un range exception --> mensaje error
	if (_layerRealID.size()>(unsigned short int)layerID) 
	{
		return _layerRealID[layerID];
	}
	else 
	{
		streamlog_out(ERROR) << "SiStripGeom::getLayerRealID - layerID: " << layerID << " out of range!!!"
			<< std::endl;
		exit(-1);
	}
}

//
// Transform real layer ID to C-type numbering 0 - n ...
//
short int SiStripGeom::getLayerIDCTypeNo(int realLayerID) const
{
std::cout << realLayerID << std::endl;
	for (unsigned int i=0; i<_layerRealID.size(); i++) 
	{
		if (realLayerID == _layerRealID[i]) 
		{
			return i;
		}
	}
	
	streamlog_out(ERROR) << "SiStripGeom::getLayerIDCTypeNo - layer: " << 
		realLayerID << " not found in layer list!!!" << std::endl;
	exit(-1);
}

//
// Get layer type
//
short int SiStripGeom::getLayerType(short int layerID) const
{
	//TODO: Posible mejora del algoritmo: try-catch
	//      capturando un range exception --> mensaje error
	if (_layerType.size()>(unsigned short int)layerID)
	{
		return _layerType[layerID];
	}
	else 
	{
		streamlog_out(ERROR) << "SiStripGeom::getLayerType - layerID: " << layerID << " out of range!!!"
			<< std::endl;
		exit(0);
	}
}

//
// Get layer radius
//
double SiStripGeom::getLayerRadius(short int layerID) const
{
	//TODO: Posible mejora del algoritmo: try-catch
	//      capturando un range exception --> mensaje error
	if (_layerRadius.size()>(unsigned short int)layerID) 
	{
		return _layerRadius[layerID];
	}
	else 
	{
		streamlog_out(ERROR) << "SiStripGeom::getLayerRadius - layerID: " << layerID << " out of range!!!"
			<< std::endl;
		exit(-1);
	}
}

//
// Get layer phi zero angle
//
double SiStripGeom::getLayerPhi0(short int layerID) const
{
	//TODO: Posible mejora del algoritmo: try-catch
	//      capturando un range exception --> mensaje error
	if (_layerPhi0.size()>(unsigned short int)layerID)
	{
		return _layerPhi0[layerID];
	}
	else 
	{
		streamlog_out(ERROR) << "SiStripGeom::getLayerPhi0 - layerID: " << layerID << " out of range!!!"
			<< std::endl;
		exit(0);
	}
}

//
// Get number of ladders
//
short int SiStripGeom::getNLadders(short int layerID) const
{
	//TODO: Posible mejora del algoritmo: try-catch
	//      capturando un range exception --> mensaje error
	if (_numberOfLadders.size()>(unsigned short int)layerID)
	{
		return _numberOfLadders[layerID];
	}
	else 
	{
		streamlog_out(ERROR) << "SiStripGeom::getNLadders - layerID: " << layerID << " out of range!!!"
                           << std::endl;
      exit(0);
   }
}

//
// Get ladder thickness
//
double SiStripGeom::getLadderThick(short int layerID) const
{
	//TODO: Posible mejora del algoritmo: try-catch
	//      capturando un range exception --> mensaje error
	if (_ladderThick.size()>(unsigned short int)layerID)
	{
		return _ladderThick[layerID];
	}
   	else 
	{
		streamlog_out(ERROR) << "SiStripGeom::getLadderThick - layerID: " << layerID << " out of range!!!"
			<< std::endl;
		exit(-1);
	}
}

// FIXME::TODO: queda pendiente desde aqui hacia abajo

//
// Get ladder width
//
double SiStripGeom::getLadderWidth(short int layerID) const
{
   if (_ladderWidth.size()>(unsigned short int)layerID) return _ladderWidth[layerID];
   else {

      streamlog_out(ERROR) << "SiStripGeom::getLadderWidth - layerID: " << layerID << " out of range!!!"
                           << std::endl;
      exit(0);
   }
}

//
// Get ladder length
//
double SiStripGeom::getLadderLength(short int layerID) const
{
   if (_ladderLength.size()>(unsigned short int)layerID) return _ladderLength[layerID];
   else {

      streamlog_out(ERROR) << "SiStripGeom::getLadderLength - layerID: " << layerID << " out of range!!!"
                           << std::endl;
      exit(0);
   }
}

//
// Get ladder offset in Y
//
double SiStripGeom::getLadderOffsetY(short int layerID) const
{
   if (_ladderOffsetY.size()>(unsigned short int)layerID) return _ladderOffsetY[layerID];
   else {

      streamlog_out(ERROR) << "SiStripGeom::getLadderOffsetY - layerID: " << layerID << " out of range!!!"
                           << std::endl;
      exit(0);
   }
}

//
// Get ladder offset in Z
//
double SiStripGeom::getLadderOffsetZ(short int layerID) const
{
   if (_ladderOffsetZ.size()>(unsigned short int)layerID) return _ladderOffsetZ[layerID];
   else {

      streamlog_out(ERROR) << "SiStripGeom::getLadderOffsetZ - layerID: " << layerID << " out of range!!!"
                           << std::endl;
      exit(0);
   }
}

//
// Get ladder rotation - phi angle (in system of units defined in PhysicalConstants.h)
//
double SiStripGeom::getLadderPhi(short int layerID, short int ladderID) const
{
   if (ladderID<getNLadders(layerID)) return (getLayerPhi0(layerID) + 2*pi/getNLadders(layerID)*ladderID);
   else {

      streamlog_out(ERROR) << "SiStripGeom::getLadderPhi - ladderID: " << ladderID << " out of range!!!"
                           << std::endl;
      exit(0);
   }
}

//
// Get ladder rotation - theta angle
//
double SiStripGeom::getLadderTheta(short int layerID) const
{
   if (_layerTheta.size()>(unsigned short int)layerID) return _layerTheta[layerID];
   else {

      streamlog_out(ERROR) << "SiStripGeom::getLadderTheta - layerID: " << layerID << " out of range!!!"
                           << std::endl;
      exit(0);
   }
}

//
// Get number of sensors for given ladder
//
short int SiStripGeom::getNSensors(short int layerID) const
{
   if (_numberOfSensors.size()>(unsigned short int)layerID) return _numberOfSensors[layerID];
   else {

      streamlog_out(ERROR) << "SiStripGeom::getNSensors - layerID: " << layerID << " out of range!!!"
                           << std::endl;
      exit(0);
   }
}

//
// Get number of strips in Z axis (in each sensor)
//
int SiStripGeom::getSensorNStripsInZ(short int layerID) const
{
   if (_sensorNStripsInZ.size()>(unsigned short int)layerID) return _sensorNStripsInZ[layerID];
   else return 0;
}

//
// Get number of strips in R-Phi (in each sensor)
//
int SiStripGeom::getSensorNStripsInRPhi(short int layerID) const
{
   if (_sensorNStripsInRPhi.size()>(unsigned short int)layerID) return _sensorNStripsInRPhi[layerID];
   else return 0;
}

//
// Get sensor pitch in Z axis for barrel-type and forward-type sensors
//
double SiStripGeom::getSensorPitchInZ(short int layerID) const
{
   if (_sensorPitchInZ.size()>(unsigned short int)layerID) return _sensorPitchInZ[layerID];
   else return 0.;
}

//
// Get sensor pitch in R-Phi for barrel-type and forward-type sensors
//
double SiStripGeom::getSensorPitchInRPhi(short int layerID, double posZ) const
{
   if      ( (getLayerType(layerID) == stripB) && (_sensorPitchInRPhi.size()>(unsigned short int)layerID) ) return _sensorPitchInRPhi[layerID];
   else if (  getLayerType(layerID) == stripF) {

      double tanAlpha = (getSensorWidth(layerID) - getSensorWidth2(layerID))/2./getSensorLength(layerID);
      if ( (posZ>-EPS*um) && (posZ<(getSensorLength(layerID) + EPS*um)) ) return ((getSensorWidth(layerID) - 2*tanAlpha*posZ)/getSensorNStripsInRPhi(layerID));
      else {

         streamlog_out(ERROR) << std::setiosflags(std::ios::fixed | std::ios::internal ) << std::setprecision(2)
                              << "SiStripGeom::getSensorPitchInRPhi - posZ: "            << std::setw(4)         << posZ << " out of range!!!"
                              << std::resetiosflags(std::ios::showpos)
                              << std::setprecision(0)
                              << std::endl;
         exit(0);;
      }
   }
   else return 0.;
}

//
// Get sensor thickness
//
double SiStripGeom::getSensorThick(short int layerID) const
{
   if (_sensorThick.size()>(unsigned short int)layerID) return _sensorThick[layerID];
   else {

      streamlog_out(ERROR) << "SiStripGeom::getSensorThick - layerID: " << layerID << " out of range!!!"
                           << std::endl;
      exit(0);
   }
}

//
// Get sensor width
//
double SiStripGeom::getSensorWidth(short int layerID) const
{
   if (_sensorWidth.size()>(unsigned short int)layerID) return _sensorWidth[layerID];
   else {

      streamlog_out(ERROR) << "SiStripGeom::getSensorWidth - layerID: " << layerID << " out of range!!!"
                           << std::endl;
      exit(0);
   }
}

//
// Get sensor width2
//
double SiStripGeom::getSensorWidth2(short int layerID) const
{
   if (_sensorWidth2.size()>(unsigned short int)layerID) return _sensorWidth2[layerID];
   else                                                  return getSensorWidth(layerID);
}

//
// Get sensor length
//
double SiStripGeom::getSensorLength(short int layerID) const
{
   if (_sensorLength.size()>(unsigned short int)layerID) return _sensorLength[layerID];
   else {

      streamlog_out(ERROR) << "SiStripGeom::getSensorLength - layerID: " << layerID << " out of range!!!"
                           << std::endl;
      exit(0);
   }
}

//
// Get gap size inbetween sensors
//
double SiStripGeom::getSensorGapInBetween(short int layerID) const
{
   if (_sensorGapInBtw.size()>(unsigned short int)layerID) return _sensorGapInBtw[layerID];
   else {

      streamlog_out(ERROR) << "SiStripGeom::getSensorGapInBetween - layerID: " << layerID << " out of range!!!"
                           << std::endl;
      exit(0);
   }
}

//
// Get width of sensor rim in Z (passive part of silicon)
//
double SiStripGeom::getSensorRimWidthInZ(short int layerID) const
{
  if (_sensorRimWidthInZ.size()>(unsigned short int)layerID) return _sensorRimWidthInZ[layerID];
   else {

      streamlog_out(ERROR) << "SiStripGeom::getSensorRimWidthInZ - layerID: " << layerID << " out of range!!!"
                           << std::endl;
      exit(0);
   } 
}

//
// Get width of sensor rim in R-Phi (passive part of silicon)
//
double SiStripGeom::getSensorRimWidthInRPhi(short int layerID) const
{
   if (_sensorRimWidthInRPhi.size()>(unsigned short int)layerID) return _sensorRimWidthInRPhi[layerID];
   else {

      streamlog_out(ERROR) << "SiStripGeom::getSensorRimWidthInRPhi - layerID: " << layerID << " out of range!!!"
                           << std::endl;
      exit(0);
   }
}


// TRANSFORMATION METHODS - GLOBAL REF. SYSTEM

//
// Method transforming given point from global ref. system to local ref. system
// (parameters: layerID, ladderID, sensorID and space point in global ref. system)
//
/*
Hep3Vector SiStripGeom::transformPointToLocal(short int layerID, short int ladderID, short int sensorID, const Hep3Vector & globalPoint)
{
   // Initialize local point
	Hep3Vector localPoint(globalPoint);

   // Gear type: VXD
   if (_gearType == "VXD") {

      // Calculate rotation angles
	   double theta = getLadderTheta(layerID);
	   double phi   = getLadderPhi(layerID, ladderID);

      // Find (0,0,0) position of local coordinate system
      Hep3Vector localOrigin(getLayerRadius(layerID), getLadderOffsetY(layerID), getLadderOffsetZ(layerID));
      localOrigin.rotateZ(+phi);

      // Perform translation - to the center of a ladder
      localPoint -= localOrigin;

      // Perform rotation - to the center of a ladder
      localPoint.rotateZ(-phi);
      localPoint.rotateY(+theta);


      // Perform translation such as X, Y, Z are positive
      localPoint += Hep3Vector(+getSensorThick(layerID)/2., +getSensorWidth(layerID)/2.,
                               +0.5*getLadderLength(layerID) - (2*sensorID + 1)*getSensorRimWidthInZ(layerID) - sensorID*getSensorGapInBetween(layerID)
                               - sensorID*getSensorLength(layerID));

      // Check if local point within sensor boundaries +- epsilon
      if (isPointOutOfSensor(layerID, localPoint)) {

         streamlog_out(ERROR) << std::setprecision(3) << "SiStripGeom::transformPointToLocal - point: "
                              << localPoint           << " is out of sensor!!!"
                              << std::setprecision(0) << std::endl;
         exit(0);
      }
   }

   // Gear type: unknown - error
   else {
      streamlog_out(ERROR) << "Unknown gear type!"
                           << std::endl;

      exit(0);
   }

   // Return space point in local ref. system
   return localPoint;
}

//
// Method transforming given vector from global ref. system to local ref. system
// (parameters: layerID, ladderID, sensorID and vector in global ref. system)
//
Hep3Vector SiStripGeom::transformVecToLocal(short int layerID, short int ladderID, const Hep3Vector & globalVec)
{
   // Initialize local vector
	Hep3Vector localVec(globalVec);

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
HepMatrix SiStripGeom::transformMatxToLocal(short int layerID, short int ladderID, const HepMatrix & globalMatrix)
{
   // Initialize local matrix 3x3 to zero values
	HepMatrix localMatrix(3,3,0);

   // Initialize rotation matrices: R, R^T (transposition)
	HepMatrix rotMatrix(3,3,0);
	HepMatrix rotMatrixT(3,3,0);

	// Gear type: VXD
	   if (_gearType == "VXD") {

	      // Calculate rotation angles
	      double theta = getLadderTheta(layerID);
	      double phi   = getLadderPhi(layerID, ladderID);

	      // Calculate rotation matrices - help
	      HepRotation rotMatrixZ(Hep3Vector(0,0,1),-phi);
	      HepRotation rotMatrixY(Hep3Vector(0,1,0),+theta);
	      HepRotation rotMatrixHelp(rotMatrixY*rotMatrixZ);

	      HepRotation rotMatrixHelpT(rotMatrixHelp.inverse());

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
Hep3Vector SiStripGeom::transformPointToGlobal(short int layerID, short int ladderID, short int sensorID, const Hep3Vector & localPoint)
{
   // Initialize global point
   Hep3Vector globalPoint(localPoint);

   // Gear type: VXD
   if (_gearType == "VXD") {

      // Calculate rotation angles
      double theta = getLadderTheta(layerID);
      double phi   = getLadderPhi(layerID, ladderID);

      // Find (0,0,0) position of local coordinate system
      Hep3Vector localOrigin(getLayerRadius(layerID), getLadderOffsetY(layerID), getLadderOffsetZ(layerID));
      localOrigin.rotateZ(+phi);

      // Perform translation - to the center of a ladder
      globalPoint -= Hep3Vector(+getSensorThick(layerID)/2., +getSensorWidth(layerID)/2.,
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
Hep3Vector SiStripGeom::transformVecToGlobal(short int layerID, short int ladderID, const Hep3Vector & localVec)
{
   // Initialize global vector
	Hep3Vector globalVec(localVec);

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
HepMatrix SiStripGeom::transformMatxToGlobal(short int layerID, short int ladderID, const HepMatrix & localMatrix)
{
   // Initialize local matrix 3x3 to zero values
   HepMatrix globalMatrix(3,3,0);

   // Initialize rotation matrices: R, R^T (transposition)
   HepMatrix rotMatrix(3,3,0);
   HepMatrix rotMatrixT(3,3,0);

   // Gear type: VXD
   if (_gearType == "VXD") {

      // Calculate rotation angles
      double theta = getLadderTheta(layerID);
      double phi   = getLadderPhi(layerID, ladderID);

      // Calculate rotation matrices - help
      HepRotation rotMatrixY(Hep3Vector(0,1,0),-theta);
      HepRotation rotMatrixZ(Hep3Vector(0,0,1),+phi);
      HepRotation rotMatrixHelp(rotMatrixZ*rotMatrixY);

      HepRotation rotMatrixHelpT(rotMatrixHelp.inverse());

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
bool SiStripGeom::isPointInsideSensor (short int layerID, short int ladderID, short int sensorID, const Hep3Vector & point) const
{
   bool isIn = false;

   // Gear type: VXD
   if (_gearType == "VXD") {
      //
      // Transform point to local ref. system

      // Initialize local point
      Hep3Vector localPoint(point);

      // Calculate rotation angles
      double theta = getLadderTheta(layerID);
      double phi   = getLadderPhi(layerID, ladderID);

      // Find (0,0,0) position of local coordinate system
      Hep3Vector localOrigin(getLayerRadius(layerID), getLadderOffsetY(layerID), getLadderOffsetZ(layerID));
      localOrigin.rotateZ(+phi);

      // Perform translation - to the center of a ladder
      localPoint -= localOrigin;

      // Perform rotation - to the center of a ladder
      localPoint.rotateZ(-phi);
      localPoint.rotateY(+theta);

      // Perform translation such as X, Y, Z are positive
      localPoint += Hep3Vector(+getSensorThick(layerID)/2., +getSensorWidth(layerID)/2.,
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
         Hep3Vector recalcPoint(localPoint.getX(), localPoint.getY() - getSensorWidth(layerID)/2. + actualSensorWidth/2., localPoint.getZ());

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
bool SiStripGeom::isPointOutOfSensor(short int layerID, const Hep3Vector & point) const
{
	bool isOut = false;

   // Gear type: VXD
   if (_gearType == "VXD") {

      // Boundary set +- epsilon
      // Barrel-type sensors
      if ( (getLayerType(layerID) == stripB) && ( (point.getX() > (getSensorThick(layerID) +EPS*um)) || (point.getX() < (-EPS*um)) ||
                                                  (point.getY() > (getSensorWidth(layerID) +EPS*um)) || (point.getY() < (-EPS*um)) ||
                                                  (point.getZ() > (getSensorLength(layerID)+EPS*um)) || (point.getZ() < (-EPS*um)) ) ) isOut = true;
      // Forward-type sensors
      if (getLayerType(layerID) == stripF) {

         double tanAlpha          = (getSensorWidth(layerID) - getSensorWidth2(layerID))/2./getSensorLength(layerID);
         double actualSensorWidth = (getSensorWidth(layerID) - 2*tanAlpha*point.getZ());

         // Recalculate point into "local" local system (width depends on posZ)
         Hep3Vector recalcPoint(point.getX(), point.getY() - getSensorWidth(layerID)/2. + actualSensorWidth/2., point.getZ());

         if ( (recalcPoint.getX() > (getSensorThick(layerID) +EPS*um)) || (recalcPoint.getX() < (-EPS*um)) ||
              (recalcPoint.getY() > (actualSensorWidth       +EPS*um)) || (recalcPoint.getY() < (-EPS*um)) ||
              (recalcPoint.getZ() > (getSensorLength(layerID)+EPS*um)) || (recalcPoint.getZ() < (-EPS*um)) ) isOut = true;
      }

   }
   // Gear type: unknown - error
   else {
   	streamlog_out(ERROR) << "Unknown gear type!"
   	                     << std::endl;

      exit(0);
   }

   // Return if out or not
   return isOut;
}

*/

//
// Get Z-position of given strip in local ref system (in system of units defined in PhysicalConstants.h);
// strips are considered to be perpendicular to beam axis for both barrel-type and forward-type sensors.
//
double SiStripGeom::getStripPosInZ(short int layerID, int stripID) const
{
	
	// Get pitch
	double sensPitch = getSensorPitchInZ(layerID);
	
	// Calculate position
	double posZ = sensPitch*(stripID + 0.5);
	// Error
	if( (posZ<0.) || (posZ>getSensorLength(layerID)) ) 
	{
		streamlog_out(ERROR) << "SiStripGeom::getStripPosInZ - position out of sensor!!!"
			<< std::endl;
		exit(-1);
	}

	// Return Z position of given strip in local ref. system
	return posZ;
}

//
// Get R-Phi position of given strip in local ref system (in system of units defined in PhysicalConstants.h);
// strips are considered to be parallel to beam axis for barrel-type sensors and at angle alpha for forward-type
// sensors (see getSensorPitchInRPhi method).
//
double SiStripGeom::getStripPosInRPhi(short int layerID, int stripID, double posZ) const
{
	// Get pitch
	double sensPitch = getSensorPitchInRPhi(layerID, posZ);
	
	// Calculate position
	double posRPhi = sensPitch*(stripID + 0.5);
	
	// Error
	if ( (posRPhi<0.) || (posRPhi>getSensorWidth(layerID)) ) 
	{
		streamlog_out(ERROR) << "SiStripGeom::getStripPosInRPhi - position out of sensor!!!"
			<< std::endl;
		exit(-1);
	}
	
	// Return R-Phi position of given strip in local ref. system
	return posRPhi;
}

//
// Get strip ID (in Z), point is given in local ref. system; strips are
// considered to be perpendicular to beam axis for both barrel-type and
// forward-type sensors.
//
int SiStripGeom::getStripIDInZ(short int layerID, double posZ ) const
{
	// Get pitch
	double sensPitch = getSensorPitchInZ(layerID);

	if (sensPitch == 0) 
	{
		streamlog_out(ERROR) << "SiStripGeom::getStripIDInZ - division by zero (sensPitch is zero)!!!"
			<< std::endl;
		exit(-1);
	}

	// Get number of strips
	int sensNStrips = getSensorNStripsInZ(layerID);

	int stripID;
	// Calculate stripID
	if (posZ <= 0.)
	{
		stripID = 0;
	}
	else 
	{
	   stripID = floor(posZ/sensPitch);
	   if (stripID >= sensNStrips)
	   {
		   stripID = sensNStrips - 1;
	   }
	}

   	// Error
	if (stripID >= sensNStrips) 
	{
		streamlog_out(ERROR) << "SiStripGeom::getStripIDInZ - stripID in Z greater than number of strips!!!"
			<< std::endl;
		exit(-1);
	}
	// Return stripID
	return stripID;
}

//
// Get strip ID (in R-Phi), point is given in local ref. system; strips are
// considered to be parallel to beam axis for barrel-type sensors and at angle
// alpha for forward-type sensors (see getSensorPitchInRPhi method).
//
int SiStripGeom::getStripIDInRPhi(short int layerID, double posRPhi, double posZ ) const
{
	// Get pitch
	double sensPitch = getSensorPitchInRPhi(layerID, posZ);
	if (sensPitch == 0) 
	{
		streamlog_out(ERROR) << "SiStripGeom::getStripIDInRPhi - division by zero (sensPitch is zero)!!!"
			<< std::endl;
		exit(-1);
	}
	// Get number of strips
	int sensNStrips = getSensorNStripsInRPhi(layerID);
	
   	int stripID;
	// Calculate stripID
	if (posRPhi <= 0.)
	{
		stripID = 0;
	}
	else 
	{
		stripID = floor(posRPhi/sensPitch);
		if (stripID >= sensNStrips) 
		{
			stripID = sensNStrips - 1;
		}
	}
	// Error
	if (stripID >= sensNStrips) 
	{
		streamlog_out(ERROR) << "SiStripGeom::getStripIDInRPhi - stripID in RPhi greater than number of strips!!!"
			<< std::endl;
		exit(-1);
	}
	// Return stripID
	return stripID;
}


// PRINT METHODS

//
// Method printing general Gear parameters
//
/*
void SiStripGeom::printGearParams() const
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
*/
//
// Method printing sensor Gear parameters (parameters: sensorID)
//
/*void SiStripGeom::printSensorParams(short int layerID) const
{

	// Gear type: VXD
   if (_gearType == "VXD") {

      // Print sensor parameters
      streamlog_out(MESSAGE2) << "    Parameters: " << _sensorThick[layerID]/um       << "um thick, "
                              << "with "            << _sensorPitchInZ[layerID]/um    << "um pitch "
                              << "and "             << _sensorNStripsInZ[layerID]     << " strips in Z"
                              << ", resp. "         << _sensorPitchInRPhi[layerID]/um << "um pitch "
                              << "and "             << _sensorNStripsInRPhi[layerID]  << " strips in R-Phi."
                              << std::endl;
   }
   // Gear type: unknown - error
   else {
   	streamlog_out(ERROR) << "Unknown gear type!"
   	                     << std::endl;

      exit(0);
   }
}
*/
} // Namespace;

