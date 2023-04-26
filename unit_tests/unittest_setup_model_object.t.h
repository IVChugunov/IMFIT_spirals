// Unit tests for code in setup_model_object.cpp

// Evidence for the idea that we're now doing integration tests is the large
// number of other modules that have to be compiled along with model_object.cpp
// to get these tests to work...

// See run_unittest_setup_model_object.sh for how to compile & run this


// We assume that the input nColumnsRowsVector has the following entries:
// nColumnsRowsVector[0] = nColumns
// nColumnsRowsVector[1] = nRows
// nColumnsRowsVector[2] = nColumns_psf  [optional]
// nColumnsRowsVector[3] = nRows_psf  [optional]
// nColumnsRowsVector[4] = nColumns_psf_oversampled  [optional]
// nColumnsRowsVector[5] = nRows_psf_oversampled  [optional]


#include <cxxtest/TestSuite.h>

#include <string>
#include <vector>
#include <memory>
#include <stdlib.h>
#include <math.h>

using namespace std;

#include "definitions.h"
#include "setup_model_object.h"
#include "model_object.h"
#include "psf_oversampling_info.h"
#include "config_file_parser.h"
#include "options_base.h"
#include "options_makeimage.h"
#include "options_imfit.h"
#include "options_mcmc.h"


#define SIMPLE_CONFIG_FILE "tests/config_imfit_flatsky.dat"
#define CONFIG_FILE "tests/config_imfit_poisson.dat"


// Reference things
const string  headerLine_correct = "# X0_1		Y0_1		PA_1	ell_1	I_0_1	h_1	I_sky_2	";

double psfPixels0[9] = {0.1, 0.1, 0.1, 0.1, 0.5, 0.1, 0.1, 0.1, 0.1};
double psfPixels1[9] = {0.1, 0.1, 0.1, 0.1, 1.0, 0.1, 0.1, 0.1, 0.1};
const int nColsPsf0 = 3;
const int nRowsPsf0 = 3;
const int nPixelsPsf0 = 9;


class NewTestSuite : public CxxTest::TestSuite 
{
public:
  int  status;
  configOptions  userConfigOptions1, userConfigOptions3;
  double  *smallDataImage;
  double  *smallErrorImage;
  double  *smallVarianceImage;
  double  *smallWeightImage;
  double  *smallMaskImage;
  double  *smallPSFImage;
  double  *oversampledPSFImage;
  int nSmallDataCols, nSmallDataRows;
  int nBigDataCols, nBigDataRows;
  int nSmallPSFCols, nSmallPSFRows;
  int nOsampPSFCols, nOsampPSFRows;


  // Note that setUp() gets called prior to *each* individual test function!
  void setUp()
  {
    int  status;
    string  filename1 = CONFIG_FILE;
    string  filename3 = SIMPLE_CONFIG_FILE;
    
    nSmallDataCols = nSmallDataRows = 2;
    nBigDataCols = nBigDataRows = 10;  // useful for oversampling cases
    nSmallPSFCols = nSmallPSFRows = 2;
    nOsampPSFCols = nOsampPSFRows = 4;
    
    smallDataImage = (double *)calloc(nSmallDataCols*nSmallDataRows, sizeof(double));
    smallDataImage[0] = 0.25;
    smallDataImage[1] = 0.25;
    smallDataImage[2] = 0.25;
    smallDataImage[3] = 1.0;

    // small image corresponding to square root of smallDataImage
    smallErrorImage = (double *)calloc(nSmallDataCols*nSmallDataRows, sizeof(double));
    smallErrorImage[0] = 0.5;
    smallErrorImage[1] = 0.5;
    smallErrorImage[2] = 0.5;
    smallErrorImage[3] = 1.0;

    smallMaskImage = (double *)calloc(nSmallDataCols*nSmallDataRows, sizeof(double));
    smallMaskImage[0] = 0.0;
    smallMaskImage[1] = 1.0;
    smallMaskImage[2] = 0.0;
    smallMaskImage[3] = 0.0;

    smallPSFImage = (double *)calloc(nSmallPSFCols*nSmallPSFRows, sizeof(double));
    smallPSFImage[0] = 0.0;
    smallPSFImage[1] = 1.0;
    smallPSFImage[2] = 0.0;
    smallPSFImage[3] = 0.0;
  }

  void tearDown()
  {
    free(smallDataImage);
    free(smallErrorImage);
    free(smallWeightImage);
    free(smallMaskImage);
    free(smallPSFImage);
  }
  
  
  void testSetupMakeimage_simple( void )
  {
    ModelObject *theModel = NULL;
    shared_ptr<OptionsBase> optionsPtr;
    vector<int> nColumnsRowsVect;
  
    optionsPtr = make_shared<MakeimageOptions>();
    nColumnsRowsVect.push_back(nSmallDataCols);
    nColumnsRowsVect.push_back(nSmallDataRows);
  
    theModel = SetupModelObject(optionsPtr, nColumnsRowsVect, NULL);
  
    long nDataVals_true = 4;
    long nDataVals = theModel->GetNDataValues();
    TS_ASSERT_EQUALS(nDataVals, nDataVals_true);
    bool psfPresent = theModel->HasPSF();
    bool oversampledPSFPresent = theModel->HasOversampledPSF();
    TS_ASSERT_EQUALS(psfPresent, false);
    TS_ASSERT_EQUALS(oversampledPSFPresent, false);

    delete theModel;
  }

  void testSetupMakeimage_withPSF( void )
  {
    ModelObject *theModel = NULL;
    shared_ptr<OptionsBase> optionsPtr;
    vector<int> nColumnsRowsVect;

    long nDataVals_true = 4;
    long nDataVals;
  
    optionsPtr = make_shared<MakeimageOptions>();
    optionsPtr->psfImagePresent = true;
    nColumnsRowsVect.push_back(nSmallDataCols);
    nColumnsRowsVect.push_back(nSmallDataRows);
    nColumnsRowsVect.push_back(nSmallPSFCols);
    nColumnsRowsVect.push_back(nSmallPSFRows);
  
    theModel = SetupModelObject(optionsPtr, nColumnsRowsVect, NULL, smallPSFImage);
  
    nDataVals = theModel->GetNDataValues();
    TS_ASSERT_EQUALS(nDataVals, nDataVals_true);
    bool psfPresent = theModel->HasPSF();
    bool oversampledPSFPresent = theModel->HasOversampledPSF();
    TS_ASSERT_EQUALS(psfPresent, true);
    TS_ASSERT_EQUALS(oversampledPSFPresent, false);

    delete theModel;
  }

  void testSetupMakeimage_withOversampledPSF( void )
  {
    ModelObject *theModel = NULL;
    shared_ptr<OptionsBase> optionsPtr;
    vector<int> nColumnsRowsVect;
    PsfOversamplingInfo *psfOsampleInfo;
    vector<PsfOversamplingInfo *> psfOsampleInfoVect;

    oversampledPSFImage = (double *)calloc(nOsampPSFCols*nOsampPSFRows, sizeof(double));
    oversampledPSFImage[0] = 0.0;
    oversampledPSFImage[1] = 1.0;
    oversampledPSFImage[2] = 0.0;
    oversampledPSFImage[3] = 0.0;

    long nDataVals_true = 100;
    long nDataVals;
  
    optionsPtr = make_shared<MakeimageOptions>();
    optionsPtr->psfImagePresent = true;
    optionsPtr->psfOversampling = true;

    // use "big" data-image dimensions to make sure PSF oversampling
    // regions can be accommodated
    nColumnsRowsVect.push_back(nBigDataCols);
    nColumnsRowsVect.push_back(nBigDataRows);
    nColumnsRowsVect.push_back(nSmallPSFCols);
    nColumnsRowsVect.push_back(nSmallPSFRows);
    
    psfOsampleInfo = new PsfOversamplingInfo(oversampledPSFImage, 2, 2, 2, "1:2,1:2");
    psfOsampleInfoVect.push_back(psfOsampleInfo);

    theModel = SetupModelObject(optionsPtr, nColumnsRowsVect, NULL, smallPSFImage,
    							NULL, NULL, psfOsampleInfoVect);
  
    nDataVals = theModel->GetNDataValues();
    TS_ASSERT_EQUALS(nDataVals, nDataVals_true);
    bool psfPresent = theModel->HasPSF();
    bool oversampledPSFPresent = theModel->HasOversampledPSF();
    TS_ASSERT_EQUALS(psfPresent, true);
    TS_ASSERT_EQUALS(oversampledPSFPresent, true);

    delete theModel;
    delete psfOsampleInfo;
  }

  void testSetupMakeimage_withMultiOversampledPSF( void )
  {
    ModelObject *theModel = NULL;
    shared_ptr<OptionsBase> optionsPtr;
    vector<int> nColumnsRowsVect;
    PsfOversamplingInfo *psfOsampleInfo1;
    PsfOversamplingInfo *psfOsampleInfo2;
    vector<PsfOversamplingInfo *> psfOsampleInfoVect;

    oversampledPSFImage = (double *)calloc(nOsampPSFCols*nOsampPSFRows, sizeof(double));
    oversampledPSFImage[0] = 0.0;
    oversampledPSFImage[1] = 1.0;
    oversampledPSFImage[2] = 0.0;
    oversampledPSFImage[3] = 0.0;

    long nDataVals_true = 100;
    long nDataVals;
  
    optionsPtr = make_shared<MakeimageOptions>();
    optionsPtr->psfImagePresent = true;
    optionsPtr->psfOversampling = true;

    // use "big" data-image dimensions to make sure PSF oversampling
    // regions can be accommodated
    nColumnsRowsVect.push_back(nBigDataCols);
    nColumnsRowsVect.push_back(nBigDataRows);
    nColumnsRowsVect.push_back(nSmallPSFCols);
    nColumnsRowsVect.push_back(nSmallPSFRows);
    
    string regionStr = "1:2,1:2";
    psfOsampleInfo1 = new PsfOversamplingInfo(oversampledPSFImage, 2, 2, 2, regionStr, 0,0, true);
    psfOsampleInfoVect.push_back(psfOsampleInfo1);
    regionStr = "3:4,3:4";
    psfOsampleInfo2 = new PsfOversamplingInfo(oversampledPSFImage, 2, 2, 2, regionStr, 0,0, false);
    psfOsampleInfoVect.push_back(psfOsampleInfo2);

    theModel = SetupModelObject(optionsPtr, nColumnsRowsVect, NULL, smallPSFImage,
    							NULL, NULL, psfOsampleInfoVect);
  
    nDataVals = theModel->GetNDataValues();
    TS_ASSERT_EQUALS(nDataVals, nDataVals_true);
    bool psfPresent = theModel->HasPSF();
    bool oversampledPSFPresent = theModel->HasOversampledPSF();
    TS_ASSERT_EQUALS(psfPresent, true);
    TS_ASSERT_EQUALS(oversampledPSFPresent, true);

    delete theModel;
    delete psfOsampleInfo1;
    delete psfOsampleInfo2;
  }

 
  void testSetupImfit_simple( void )
  {
    ModelObject *theModel = NULL;
    shared_ptr<OptionsBase> optionsPtr;
    vector<int> nColumnsRowsVect;
  
    optionsPtr = make_shared<ImfitOptions>();
    nColumnsRowsVect.push_back(nSmallDataCols);
    nColumnsRowsVect.push_back(nSmallDataRows);
  
    theModel = SetupModelObject(optionsPtr, nColumnsRowsVect, smallDataImage);
  
    long nDataVals_true = 4;
    long nDataVals = theModel->GetNDataValues();
    TS_ASSERT_EQUALS(nDataVals, nDataVals_true);
    bool psfPresent = theModel->HasPSF();
    bool oversampledPSFPresent = theModel->HasOversampledPSF();
    TS_ASSERT_EQUALS(psfPresent, false);
    TS_ASSERT_EQUALS(oversampledPSFPresent, false);

    delete theModel;
  }

  void testSetupImfit_withPSF( void )
  {
    ModelObject *theModel = NULL;
    shared_ptr<OptionsBase> optionsPtr;
    vector<int> nColumnsRowsVect;
  
    optionsPtr = make_shared<ImfitOptions>();
    optionsPtr->psfImagePresent = true;
    nColumnsRowsVect.push_back(nSmallDataCols);
    nColumnsRowsVect.push_back(nSmallDataRows);
    nColumnsRowsVect.push_back(nSmallPSFCols);
    nColumnsRowsVect.push_back(nSmallPSFRows);
  
    theModel = SetupModelObject(optionsPtr, nColumnsRowsVect, smallDataImage,
    			smallPSFImage);
  
    long nDataVals_true = 4;
    long nDataVals = theModel->GetNDataValues();
    TS_ASSERT_EQUALS(nDataVals, nDataVals_true);
    bool psfPresent = theModel->HasPSF();
    bool oversampledPSFPresent = theModel->HasOversampledPSF();
    bool maskPresent = theModel->HasMask();
    TS_ASSERT_EQUALS(psfPresent, true);
    TS_ASSERT_EQUALS(oversampledPSFPresent, false);
    TS_ASSERT_EQUALS(maskPresent, false);

    delete theModel;
  }

  void testSetupImfit_withOversampledPSF( void )
  {
    ModelObject *theModel = NULL;
    shared_ptr<OptionsBase> optionsPtr;
    vector<int> nColumnsRowsVect;
    PsfOversamplingInfo *psfOsampleInfo;
    vector<PsfOversamplingInfo *> psfOsampleInfoVect;

    oversampledPSFImage = (double *)calloc(nOsampPSFCols*nOsampPSFRows, sizeof(double));
    oversampledPSFImage[0] = 0.0;
    oversampledPSFImage[1] = 1.0;
    oversampledPSFImage[2] = 0.0;
    oversampledPSFImage[3] = 0.0;

    optionsPtr = make_shared<ImfitOptions>();
    optionsPtr->psfImagePresent = true;
    optionsPtr->psfOversampling = true;

    nColumnsRowsVect.push_back(nSmallDataCols);
    nColumnsRowsVect.push_back(nSmallDataRows);
    nColumnsRowsVect.push_back(nSmallPSFCols);
    nColumnsRowsVect.push_back(nSmallPSFRows);

    psfOsampleInfo = new PsfOversamplingInfo(oversampledPSFImage, 2, 2, 2, "1:2,1:2");
    psfOsampleInfoVect.push_back(psfOsampleInfo);

    theModel = SetupModelObject(optionsPtr, nColumnsRowsVect, smallDataImage, smallPSFImage,
    					NULL, NULL, psfOsampleInfoVect);
  
    long nDataVals_true = 4;
    long nDataVals = theModel->GetNDataValues();
    TS_ASSERT_EQUALS(nDataVals, nDataVals_true);
    bool psfPresent = theModel->HasPSF();
    bool oversampledPSFPresent = theModel->HasOversampledPSF();
    bool maskPresent = theModel->HasMask();
    TS_ASSERT_EQUALS(psfPresent, true);
    TS_ASSERT_EQUALS(oversampledPSFPresent, true);
    TS_ASSERT_EQUALS(maskPresent, false);

    delete theModel;
    delete psfOsampleInfo;
  }

  void testSetupImfit_withMultipleOversampledPSF( void )
  {
    ModelObject *theModel = NULL;
    shared_ptr<OptionsBase> optionsPtr;
    vector<int> nColumnsRowsVect;
    PsfOversamplingInfo *psfOsampleInfo1;
    PsfOversamplingInfo *psfOsampleInfo2;
    vector<PsfOversamplingInfo *> psfOsampleInfoVect;

    double *bigDataImage = (double *)calloc(nBigDataCols*nBigDataRows, sizeof(double));
    
    oversampledPSFImage = (double *)calloc(nOsampPSFCols*nOsampPSFRows, sizeof(double));
    oversampledPSFImage[0] = 0.0;
    oversampledPSFImage[1] = 1.0;
    oversampledPSFImage[2] = 0.0;
    oversampledPSFImage[3] = 0.0;

    optionsPtr = make_shared<ImfitOptions>();
    optionsPtr->psfImagePresent = true;
    optionsPtr->psfOversampling = true;

    // use "big" data-image dimensions to make sure PSF oversampling
    // regions can be accommodated
    nColumnsRowsVect.push_back(nBigDataCols);
    nColumnsRowsVect.push_back(nBigDataRows);
    nColumnsRowsVect.push_back(nSmallPSFCols);
    nColumnsRowsVect.push_back(nSmallPSFRows);

    string regionStr = "1:2,1:2";
    psfOsampleInfo1 = new PsfOversamplingInfo(oversampledPSFImage, 2, 2, 2, regionStr, 0,0, true);
    psfOsampleInfoVect.push_back(psfOsampleInfo1);
    regionStr = "3:4,3:4";
    psfOsampleInfo2 = new PsfOversamplingInfo(oversampledPSFImage, 2, 2, 2, regionStr, 0,0, false);
    psfOsampleInfoVect.push_back(psfOsampleInfo2);

    theModel = SetupModelObject(optionsPtr, nColumnsRowsVect, bigDataImage, smallPSFImage,
    					NULL, NULL, psfOsampleInfoVect);
  
    long nDataVals_true = 100;
    long nDataVals = theModel->GetNDataValues();
    TS_ASSERT_EQUALS(nDataVals, nDataVals_true);
    bool psfPresent = theModel->HasPSF();
    bool oversampledPSFPresent = theModel->HasOversampledPSF();
    bool maskPresent = theModel->HasMask();
    TS_ASSERT_EQUALS(psfPresent, true);
    TS_ASSERT_EQUALS(oversampledPSFPresent, true);
    TS_ASSERT_EQUALS(maskPresent, false);

    delete theModel;
    delete psfOsampleInfo1;
    delete psfOsampleInfo2;
    free(bigDataImage);
  }

  void testSetupImfit_withMask( void )
  {
    ModelObject *theModel = NULL;
    shared_ptr<OptionsBase> optionsPtr;
    vector<int> nColumnsRowsVect;
  
    optionsPtr = make_shared<ImfitOptions>();
    optionsPtr->maskImagePresent = true;
    nColumnsRowsVect.push_back(nSmallDataCols);
    nColumnsRowsVect.push_back(nSmallDataRows);
  
    theModel = SetupModelObject(optionsPtr, nColumnsRowsVect, smallDataImage,
    			NULL, smallMaskImage);
  
    long nDataVals_true = 4;
    long nDataVals = theModel->GetNDataValues();
    TS_ASSERT_EQUALS(nDataVals, nDataVals_true);
    bool psfPresent = theModel->HasPSF();
    bool oversampledPSFPresent = theModel->HasOversampledPSF();
    bool maskPresent = theModel->HasMask();
    TS_ASSERT_EQUALS(psfPresent, false);
    TS_ASSERT_EQUALS(oversampledPSFPresent, false);
    TS_ASSERT_EQUALS(maskPresent, true);

    delete theModel;
  }

  void testSetupImfit_withError( void )
  {
    ModelObject *theModel = NULL;
    shared_ptr<OptionsBase> optionsPtr;
    vector<int> nColumnsRowsVect;
  
    optionsPtr = make_shared<ImfitOptions>();
    optionsPtr->noiseImagePresent = true;
    nColumnsRowsVect.push_back(nSmallDataCols);
    nColumnsRowsVect.push_back(nSmallDataRows);
  
    theModel = SetupModelObject(optionsPtr, nColumnsRowsVect, smallDataImage,
    			NULL, NULL, smallErrorImage);
  
    long nDataVals_true = 4;
    long nDataVals = theModel->GetNDataValues();
    TS_ASSERT_EQUALS(nDataVals, nDataVals_true);
    bool psfPresent = theModel->HasPSF();
    bool oversampledPSFPresent = theModel->HasOversampledPSF();
    bool maskPresent = theModel->HasMask();
    TS_ASSERT_EQUALS(psfPresent, false);
    TS_ASSERT_EQUALS(oversampledPSFPresent, false);
    TS_ASSERT_EQUALS(maskPresent, false);

    delete theModel;
  }


  void testSetupMCMC_simple( void )
  {
    ModelObject *theModel = NULL;
    shared_ptr<OptionsBase> optionsPtr;
    vector<int> nColumnsRowsVect;
  
    optionsPtr = make_shared<MCMCOptions>();
    nColumnsRowsVect.push_back(nSmallDataCols);
    nColumnsRowsVect.push_back(nSmallDataRows);
  
    theModel = SetupModelObject(optionsPtr, nColumnsRowsVect, smallDataImage);
  
    long nDataVals_true = 4;
    long nDataVals = theModel->GetNDataValues();
    TS_ASSERT_EQUALS(nDataVals, nDataVals_true);
    bool psfPresent = theModel->HasPSF();
    bool oversampledPSFPresent = theModel->HasOversampledPSF();
    bool maskPresent = theModel->HasMask();
    TS_ASSERT_EQUALS(psfPresent, false);
    TS_ASSERT_EQUALS(oversampledPSFPresent, false);
    TS_ASSERT_EQUALS(maskPresent, false);
	
    delete theModel;
  }

  void testSetupMCMC_withPSF( void )
  {
    ModelObject *theModel = NULL;
    shared_ptr<OptionsBase> optionsPtr;
    vector<int> nColumnsRowsVect;
  
    optionsPtr = make_shared<MCMCOptions>();
    optionsPtr->psfImagePresent = true;
    nColumnsRowsVect.push_back(nSmallDataCols);
    nColumnsRowsVect.push_back(nSmallDataRows);
    nColumnsRowsVect.push_back(nSmallPSFCols);
    nColumnsRowsVect.push_back(nSmallPSFRows);
  
    theModel = SetupModelObject(optionsPtr, nColumnsRowsVect, smallDataImage,
    			smallPSFImage);
  
    long nDataVals_true = 4;
    long nDataVals = theModel->GetNDataValues();
    TS_ASSERT_EQUALS(nDataVals, nDataVals_true);
    bool psfPresent = theModel->HasPSF();
    bool oversampledPSFPresent = theModel->HasOversampledPSF();
    bool maskPresent = theModel->HasMask();
    TS_ASSERT_EQUALS(psfPresent, true);
    TS_ASSERT_EQUALS(oversampledPSFPresent, false);
    TS_ASSERT_EQUALS(maskPresent, false);

    delete theModel;
  }

  void testSetupMCMC_withOversampledPSF( void )
  {
    ModelObject *theModel = NULL;
    shared_ptr<OptionsBase> optionsPtr;
    vector<int> nColumnsRowsVect;
    PsfOversamplingInfo *psfOsampleInfo1;
    PsfOversamplingInfo *psfOsampleInfo2;
    vector<PsfOversamplingInfo *> psfOsampleInfoVect;

    double *bigDataImage = (double *)calloc(nBigDataCols*nBigDataRows, sizeof(double));
    
    oversampledPSFImage = (double *)calloc(nOsampPSFCols*nOsampPSFRows, sizeof(double));
    oversampledPSFImage[0] = 0.0;
    oversampledPSFImage[1] = 1.0;
    oversampledPSFImage[2] = 0.0;
    oversampledPSFImage[3] = 0.0;

    optionsPtr = make_shared<MCMCOptions>();
    optionsPtr->psfImagePresent = true;
    optionsPtr->psfOversampling = true;

    nColumnsRowsVect.push_back(nBigDataCols);
    nColumnsRowsVect.push_back(nBigDataRows);
    nColumnsRowsVect.push_back(nSmallPSFCols);
    nColumnsRowsVect.push_back(nSmallPSFRows);
    
    string regionStr = "1:2,1:2";
    psfOsampleInfo1 = new PsfOversamplingInfo(oversampledPSFImage, 2, 2, 2, regionStr, 0,0, true);
    psfOsampleInfoVect.push_back(psfOsampleInfo1);
    regionStr = "3:4,3:4";
    psfOsampleInfo2 = new PsfOversamplingInfo(oversampledPSFImage, 2, 2, 2, regionStr, 0,0, false);
    psfOsampleInfoVect.push_back(psfOsampleInfo2);
  
    theModel = SetupModelObject(optionsPtr, nColumnsRowsVect, smallDataImage, smallPSFImage,
		    					NULL, NULL, psfOsampleInfoVect);
  
    long nDataVals_true = 100;
    long nDataVals = theModel->GetNDataValues();
    TS_ASSERT_EQUALS(nDataVals, nDataVals_true);
    bool psfPresent = theModel->HasPSF();
    bool oversampledPSFPresent = theModel->HasOversampledPSF();
    bool maskPresent = theModel->HasMask();
    TS_ASSERT_EQUALS(psfPresent, true);
    TS_ASSERT_EQUALS(oversampledPSFPresent, true);
    TS_ASSERT_EQUALS(maskPresent, false);

    delete theModel;
    delete psfOsampleInfo1;
    delete psfOsampleInfo2;
    free(bigDataImage);
  }

  void testSetupMCMC_withMask( void )
  {
    ModelObject *theModel = NULL;
    shared_ptr<OptionsBase> optionsPtr;
    vector<int> nColumnsRowsVect;
  
    optionsPtr = make_shared<MCMCOptions>();
    optionsPtr->maskImagePresent = true;
    nColumnsRowsVect.push_back(nSmallDataCols);
    nColumnsRowsVect.push_back(nSmallDataRows);
  
    theModel = SetupModelObject(optionsPtr, nColumnsRowsVect, smallDataImage,
    			NULL, smallMaskImage);
  
    long nDataVals_true = 4;
    long nDataVals = theModel->GetNDataValues();
    TS_ASSERT_EQUALS(nDataVals, nDataVals_true);
    bool psfPresent = theModel->HasPSF();
    bool oversampledPSFPresent = theModel->HasOversampledPSF();
    bool maskPresent = theModel->HasMask();
    TS_ASSERT_EQUALS(psfPresent, false);
    TS_ASSERT_EQUALS(oversampledPSFPresent, false);
    TS_ASSERT_EQUALS(maskPresent, true);

    delete theModel;
  }

  void testSetupMCMC_withError( void )
  {
    ModelObject *theModel = NULL;
    shared_ptr<OptionsBase> optionsPtr;
    vector<int> nColumnsRowsVect;
  
    optionsPtr = make_shared<MCMCOptions>();
    optionsPtr->noiseImagePresent = true;
    nColumnsRowsVect.push_back(nSmallDataCols);
    nColumnsRowsVect.push_back(nSmallDataRows);
  
    theModel = SetupModelObject(optionsPtr, nColumnsRowsVect, smallDataImage,
    			NULL, NULL, smallErrorImage);
  
    long nDataVals_true = 4;
    long nDataVals = theModel->GetNDataValues();
    TS_ASSERT_EQUALS(nDataVals, nDataVals_true);
    bool psfPresent = theModel->HasPSF();
    bool oversampledPSFPresent = theModel->HasOversampledPSF();
    bool maskPresent = theModel->HasMask();
    TS_ASSERT_EQUALS(psfPresent, false);
    TS_ASSERT_EQUALS(oversampledPSFPresent, false);
    TS_ASSERT_EQUALS(maskPresent, false);

    delete theModel;
  }

};
