#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkImageDuplicator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkNeighborhoodIterator.h"

const int MAX_VAL=100000;

typedef short datatype;
typedef itk::Image<datatype, 3> ImageType;
typedef itk::ImageFileReader<ImageType> ReaderType;
typedef itk::ImageFileWriter<ImageType> WriterType;
typedef itk::ImageRegionIteratorWithIndex<ImageType> IteratorType;
//static void CreateImage(ImageType::Pointer image);
 
//Expects smart pointers to scalar image objects. These need to be pointing to 4 separate memory blocks
void geos(ImageType::Pointer Inp, ImageType::Pointer Init, ImageType::Pointer Geos, ImageType::Pointer Gamma)
{ 
	typedef itk::ImageRegionIteratorWithIndex <ImageType> IteratorType;
	IteratorType InpIt(Inp,Inp->GetRequestedRegion());
	IteratorType InitIt(Init,Init->GetRequestedRegion());
	IteratorType GeosIt(Geos,Geos->GetRequestedRegion());
	IteratorType GamIt(Gamma,Gamma->GetRequestedRegion());
	//Geodesic code goes here Initialize gamma to fixed value - can be changed in the future to incorporate priors
	GamIt.GoToBegin();
	while(!GamIt.IsAtEnd())
	{
		GamIt.Set(1.0);
		++GamIt;
	}
	//Set initial infinities
	GeosIt.GoToBegin();
	InitIt.GoToBegin();
	while(!GeosIt.IsAtEnd())
	{
		  if (InitIt.Get()!=0){
          GeosIt.Set(MAX_VAL);
		  }
		
		++GeosIt;
		++InitIt;
	}
	//Setting up the neighborhood iterator
	 ImageType::SizeType radius;
     radius[0] = 1;
     radius[1] = 1;
     radius[2] = 1;
	 itk::NeighborhoodIterator<ImageType> ResNIt(radius, Geos, Geos->GetRequestedRegion());
	 itk::NeighborhoodIterator<ImageType> InpNIt(radius, Inp, Inp->GetRequestedRegion());

	//The main loops
	std::cout<<"Main loops execution: Forward pass"<<std::endl;
	ResNIt.GoToBegin();
	InpNIt.GoToBegin();
	GamIt.GoToBegin();
	InpIt.GoToBegin();
	GeosIt.GoToBegin();
	InitIt.GoToBegin();
	while(!InpIt.IsAtEnd())
	{ 
	  double C_f_arr[14];
      C_f_arr[13]=GeosIt.Get();
      C_f_arr[0]=ResNIt.GetPixel(4)+sqrt(1.0+GamIt.Get()*((InpNIt.GetPixel(13)-InpNIt.GetPixel(4))*(InpNIt.GetPixel(13)-InpNIt.GetPixel(4))));
      C_f_arr[1]=ResNIt.GetPixel(10)+sqrt(1.0+GamIt.Get()*((InpNIt.GetPixel(13)-InpNIt.GetPixel(10))*(InpNIt.GetPixel(13)-InpNIt.GetPixel(10))));
      C_f_arr[2]=ResNIt.GetPixel(12)+sqrt(1.0+GamIt.Get()*((InpNIt.GetPixel(13)-InpNIt.GetPixel(12))*(InpNIt.GetPixel(13)-InpNIt.GetPixel(12))));
      C_f_arr[3]=ResNIt.GetPixel(1)+sqrt(2.0+GamIt.Get()*((InpNIt.GetPixel(13)-InpNIt.GetPixel(1))*(InpNIt.GetPixel(13)-InpNIt.GetPixel(1))));
      C_f_arr[4]=ResNIt.GetPixel(3)+sqrt(2.0+GamIt.Get()*((InpNIt.GetPixel(13)-InpNIt.GetPixel(3))*(InpNIt.GetPixel(13)-InpNIt.GetPixel(3))));
      C_f_arr[5]=ResNIt.GetPixel(9)+sqrt(2.0+GamIt.Get()*((InpNIt.GetPixel(13)-InpNIt.GetPixel(9))*(InpNIt.GetPixel(13)-InpNIt.GetPixel(9))));
      C_f_arr[6]=ResNIt.GetPixel(0)+sqrt(3.0+GamIt.Get()*((InpNIt.GetPixel(13)-InpNIt.GetPixel(0))*(InpNIt.GetPixel(13)-InpNIt.GetPixel(0))));
      C_f_arr[7]=ResNIt.GetPixel(7)+sqrt(2.0+GamIt.Get()*((InpNIt.GetPixel(13)-InpNIt.GetPixel(7))*(InpNIt.GetPixel(13)-InpNIt.GetPixel(7))));
      C_f_arr[8]=ResNIt.GetPixel(6)+sqrt(3.0+GamIt.Get()*((InpNIt.GetPixel(13)-InpNIt.GetPixel(6))*(InpNIt.GetPixel(13)-InpNIt.GetPixel(6))));
      C_f_arr[9]=ResNIt.GetPixel(15)+sqrt(2.0+GamIt.Get()*((InpNIt.GetPixel(13)-InpNIt.GetPixel(15))*(InpNIt.GetPixel(13)-InpNIt.GetPixel(15))));
      C_f_arr[10]=ResNIt.GetPixel(24)+sqrt(3.0+GamIt.Get()*((InpNIt.GetPixel(13)-InpNIt.GetPixel(24))*(InpNIt.GetPixel(13)-InpNIt.GetPixel(24))));
      C_f_arr[11]=ResNIt.GetPixel(21)+sqrt(2.0+GamIt.Get()*((InpNIt.GetPixel(13)-InpNIt.GetPixel(21))*(InpNIt.GetPixel(13)-InpNIt.GetPixel(21))));
      C_f_arr[12]=ResNIt.GetPixel(18)+sqrt(3.0+GamIt.Get()*((InpNIt.GetPixel(13)-InpNIt.GetPixel(18))*(InpNIt.GetPixel(13)-InpNIt.GetPixel(18))));
      //std::cout<<InpNIt.GetPixel(13);
	  double minval=20000000;
      for (int i=0; i<14; ++i)
	  {
        if (C_f_arr[i]<minval)
          minval=C_f_arr[i];
      }
	  //if(minval>0)
	  //{std::cout<<minval<<std::endl;}
	  //std::cout<<GeosIt.GetIndex()<<std::endl;
      GeosIt.Set(minval);
	  ++InpIt;
      ++GeosIt;
      ++ResNIt;
      ++InpNIt;
      ++GamIt;
   	}
    InpIt.GoToReverseBegin();
    GeosIt.GoToReverseBegin();
    GamIt.GoToReverseBegin();
    ResNIt.GoToEnd();
    --ResNIt;
    InpNIt.GoToEnd();
    --InpNIt;
	std::cout<<"Main loops execution: Backward pass"<<std::endl;
	while(!InpIt.IsAtReverseEnd())
    {
	   double C_b_arr[14];
       C_b_arr[13]=GeosIt.Get();
       C_b_arr[0]=ResNIt.GetPixel(22)+sqrt(1.0+GamIt.Get()*((InpNIt.GetPixel(13)-InpNIt.GetPixel(22))*(InpNIt.GetPixel(13)-InpNIt.GetPixel(22))));
       C_b_arr[1]=ResNIt.GetPixel(16)+sqrt(1.0+GamIt.Get()*((InpNIt.GetPixel(13)-InpNIt.GetPixel(16))*(InpNIt.GetPixel(13)-InpNIt.GetPixel(16))));
       C_b_arr[2]=ResNIt.GetPixel(14)+sqrt(1.0+GamIt.Get()*((InpNIt.GetPixel(13)-InpNIt.GetPixel(14))*(InpNIt.GetPixel(13)-InpNIt.GetPixel(14))));
       C_b_arr[3]=ResNIt.GetPixel(25)+sqrt(2.0+GamIt.Get()*((InpNIt.GetPixel(13)-InpNIt.GetPixel(25))*(InpNIt.GetPixel(13)-InpNIt.GetPixel(25))));
       C_b_arr[4]=ResNIt.GetPixel(23)+sqrt(2.0+GamIt.Get()*((InpNIt.GetPixel(13)-InpNIt.GetPixel(23))*(InpNIt.GetPixel(13)-InpNIt.GetPixel(23))));
       C_b_arr[5]=ResNIt.GetPixel(17)+sqrt(2.0+GamIt.Get()*((InpNIt.GetPixel(13)-InpNIt.GetPixel(17))*(InpNIt.GetPixel(13)-InpNIt.GetPixel(17))));
       C_b_arr[6]=ResNIt.GetPixel(26)+sqrt(3.0+GamIt.Get()*((InpNIt.GetPixel(13)-InpNIt.GetPixel(26))*(InpNIt.GetPixel(13)-InpNIt.GetPixel(26))));
       C_b_arr[7]=ResNIt.GetPixel(19)+sqrt(2.0+GamIt.Get()*((InpNIt.GetPixel(13)-InpNIt.GetPixel(19))*(InpNIt.GetPixel(13)-InpNIt.GetPixel(19))));
       C_b_arr[8]=ResNIt.GetPixel(20)+sqrt(3.0+GamIt.Get()*((InpNIt.GetPixel(13)-InpNIt.GetPixel(20))*(InpNIt.GetPixel(13)-InpNIt.GetPixel(20))));
       C_b_arr[9]=ResNIt.GetPixel(11)+sqrt(2.0+GamIt.Get()*((InpNIt.GetPixel(13)-InpNIt.GetPixel(11))*(InpNIt.GetPixel(13)-InpNIt.GetPixel(11))));
       C_b_arr[10]=ResNIt.GetPixel(2)+sqrt(3.0+GamIt.Get()*((InpNIt.GetPixel(13)-InpNIt.GetPixel(2))*(InpNIt.GetPixel(13)-InpNIt.GetPixel(2))));
       C_b_arr[11]=ResNIt.GetPixel(5)+sqrt(2.0+GamIt.Get()*((InpNIt.GetPixel(13)-InpNIt.GetPixel(5))*(InpNIt.GetPixel(13)-InpNIt.GetPixel(5))));
       C_b_arr[12]=ResNIt.GetPixel(8)+sqrt(3.0+GamIt.Get()*((InpNIt.GetPixel(13)-InpNIt.GetPixel(8))*(InpNIt.GetPixel(13)-InpNIt.GetPixel(8))));
       double minval=20000000;
       for (int i=0; i<14; ++i)
	   {
         if (C_b_arr[i]<minval)
           minval=C_b_arr[i];
       }
	   //std::cout<<minval;
       GeosIt.Set(minval);
	  --InpIt;
      --GeosIt;
      --ResNIt;
      --InpNIt;
      --GamIt;
   }
	for(GeosIt.GoToBegin();!GeosIt.IsAtEnd();++GeosIt)
	{
		if(GeosIt.Get()>=MAX_VAL)
		{
         GeosIt.Set(0);
        }
	}

}

int main(int argc, char *argv[])
{ //TODO: Implement error check
  // Get input image
  ReaderType::Pointer input = ReaderType::New();
  //input->SetFileName("C:\Users\bilwaj\Desktop\GeodesicFilter\Debug\1.nii");
  input->SetFileName(argv[1]);
  input->Update();
  std::cout << "1st image read.\n";
  //Assign Inp pointer to input image in memory
  ImageType::Pointer Inp=ImageType::New();
  Inp=input->GetOutput();

  //Get Initialization 
  ReaderType::Pointer initial =ReaderType::New();
  initial->SetFileName(argv[2]);
  //initial->SetFileName("C:\Users\bilwaj\Desktop\GeodesicFilter\Debug\init_spine.nii");
  initial->Update();
  std::cout << "2nd image read.\n";
  //Assign init pointer to initialization
  ImageType::Pointer Init = ImageType::New();
  Init=initial->GetOutput();

  //Copy information if you need to keep the original image unharmed
  //Alternatively the filter may be run on the original image directly
  ImageType::Pointer Geos = ImageType::New();;
  Geos->CopyInformation(Inp);
  Geos->SetRequestedRegion( Inp->GetRequestedRegion() );
  Geos->SetBufferedRegion( Inp->GetBufferedRegion() );
  Geos->Allocate();
  Geos->FillBuffer(0);

  //For future use to specify prior probabilities in segmentation as done in the SPIE paper
  ImageType::Pointer Gamma = ImageType::New();;
  Gamma->CopyInformation(Inp);
  Gamma->SetRequestedRegion( Inp->GetRequestedRegion() );
  Gamma->SetBufferedRegion( Inp->GetBufferedRegion() );
  Gamma->Allocate();
  Gamma->FillBuffer(0);

  //Speed up calculations by focusing algorithm on a smaller region than the entire image - Here non zero pixels are being used 
  ImageType::Pointer Mask = ImageType::New();;
  Mask->CopyInformation(Inp);
  Mask->SetRequestedRegion( Inp->GetRequestedRegion() );
  Mask->SetBufferedRegion( Inp->GetBufferedRegion() );
  Mask->Allocate();
  Mask->FillBuffer(0);

  //Calling the geodesic function
  geos(Inp,Init,Geos,Gamma);

  std::cout<<"Writing output"<<std::endl;


  WriterType::Pointer output=WriterType::New();
  output->SetFileName("geos.nii");
  output->SetInput(Geos);
  output->Update();

  //std::cout<<max_val<<std::endl;
  return EXIT_SUCCESS;
}
