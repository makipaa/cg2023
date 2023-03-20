///////////////////////////////////////////////////////////////////////////////
//
//      TargaImage.cpp                          Author:     Stephen Chenney
//                                              Modified:   Eric McDaniel
//                                              Date:       Fall 2004
//
//      Implementation of TargaImage methods.  You must implement the image
//  modification functions.
//
///////////////////////////////////////////////////////////////////////////////

#include "Globals.h"
#include "TargaImage.h"
#include "libtarga.h"
#include <stdlib.h>
#include <assert.h>
#include <memory.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <map>
#include <tuple>


using namespace std;

// constants
const int           RED             = 0;                // red channel
const int           GREEN           = 1;                // green channel
const int           BLUE            = 2;                // blue channel
const unsigned char BACKGROUND[3]   = { 0, 0, 0 };      // background color


// Computes n choose s, efficiently
double Binomial(int n, int s)
{
    double        res;

    res = 1;
    for (int i = 1 ; i <= s ; i++)
        res = (n - i + 1) * res / i ;

    return res;
}// Binomial


///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage() : width(0), height(0), data(NULL)
{}// TargaImage

///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(int w, int h) : width(w), height(h)
{
   data = new unsigned char[width * height * 4];
   ClearToBlack();
}// TargaImage



///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables to values given.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(int w, int h, unsigned char *d)
{
    int i;

    width = w;
    height = h;
    data = new unsigned char[width * height * 4];

    for (i = 0; i < width * height * 4; i++)
	    data[i] = d[i];
}// TargaImage

///////////////////////////////////////////////////////////////////////////////
//
//      Copy Constructor.  Initialize member to that of input
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(const TargaImage& image) 
{
   width = image.width;
   height = image.height;
   data = NULL; 
   if (image.data != NULL) {
      data = new unsigned char[width * height * 4];
      memcpy(data, image.data, sizeof(unsigned char) * width * height * 4);
   }
}


///////////////////////////////////////////////////////////////////////////////
//
//      Destructor.  Free image memory.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::~TargaImage()
{
    if (data)
        delete[] data;
}// ~TargaImage


///////////////////////////////////////////////////////////////////////////////
//
//      Converts an image to RGB form, and returns the rgb pixel data - 24 
//  bits per pixel. The returned space should be deleted when no longer 
//  required.
//
///////////////////////////////////////////////////////////////////////////////
unsigned char* TargaImage::To_RGB(void)
{
    unsigned char   *rgb = new unsigned char[width * height * 3];
    int		    i, j;

    if (! data)
	    return NULL;

    // Divide out the alpha
    for (i = 0 ; i < height ; i++)
    {
	    int in_offset = i * width * 4;
	    int out_offset = i * width * 3;

	    for (j = 0 ; j < width ; j++)
        {
	        RGBA_To_RGB(data + (in_offset + j*4), rgb + (out_offset + j*3));
	    }
    }

    return rgb;
}// TargaImage


///////////////////////////////////////////////////////////////////////////////
//
//      Save the image to a targa file. Returns 1 on success, 0 on failure.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Save_Image(const char *filename)
{
    TargaImage	*out_image = Reverse_Rows();

    if (! out_image)
	    return false;

    if (!tga_write_raw(filename, width, height, out_image->data, TGA_TRUECOLOR_32))
    {
	    cout << "TGA Save Error: %s\n", tga_error_string(tga_get_last_error());
	    return false;
    }

    delete out_image;

    return true;
}// Save_Image


///////////////////////////////////////////////////////////////////////////////
//
//      Load a targa image from a file.  Return a new TargaImage object which 
//  must be deleted by caller.  Return NULL on failure.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage* TargaImage::Load_Image(char *filename)
{
    unsigned char   *temp_data;
    TargaImage	    *temp_image;
    TargaImage	    *result;
    int		        width, height;

    if (!filename)
    {
        cout << "No filename given." << endl;
        return NULL;
    }// if

    temp_data = (unsigned char*)tga_load(filename, &width, &height, TGA_TRUECOLOR_32);
    if (!temp_data)
    {
        cout << "TGA Error: %s\n", tga_error_string(tga_get_last_error());
	    width = height = 0;
	    return NULL;
    }
    temp_image = new TargaImage(width, height, temp_data);
    free(temp_data);

    result = temp_image->Reverse_Rows();

    delete temp_image;

    return result;
}// Load_Image


///////////////////////////////////////////////////////////////////////////////
//
//      Convert image to grayscale.  Red, green, and blue channels should all 
//  contain grayscale value.  Alpha channel shoould be left unchanged.  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::To_Grayscale()
{
  
    for(int i = 0; i < width * height * 4; i += 4)
    {
        int y = (int) (data[i + RED]*0.299 + data[i + GREEN]*0.587 + data[i + BLUE]*0.114);
            
        data[i + RED] = y;
        data[i + GREEN] = y;
        data[i + BLUE] = y;
    }

	return true;
}// To_Grayscale


///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using uniform quantization.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Quant_Uniform()
{
   for(int i = 0; i < width * height * 4; i += 4)
   {
        data[i + RED] = (int) floor(data[i + RED] / 32)*32;
        data[i + GREEN] = (int) floor(data[i + GREEN] / 32)*32;
        data[i + BLUE] = (int) floor(data[i + BLUE] / 64)*64;
   }

    return true;
}// Quant_Uniform

bool compare_int_desc(const pair<tuple<int, int, int>, int> i, const pair<tuple<int, int, int>, int> j)
{
    return (i.second > j.second);
}


vector<pair<tuple<int, int, int>, int>> get_first_256_colors_by_popularity(int height, int width, unsigned char *data)
{
    map<tuple<int, int, int>, int> hist = {};

    for(int i = 0; i < width * height * 4; i += 4)
    { 
        int r = data[i + RED];
        int g = data[i + GREEN];
        int b = data[i + BLUE];

        tuple<int, int, int> rgb = make_tuple(r, g, b);

        map<tuple<int, int, int>, int>::iterator it = hist.find(rgb);

        // Value found, increment the count
        if(it != hist.end()) {
            it->second++;
        } else {
            hist.insert(std::make_pair(rgb, 1));
        }
    }

    // Insert the map values to vector
    vector<pair<tuple<int, int, int>, int>> sorted_hist;
    map<tuple<int, int, int>, int>::iterator it2;

    for(it2 = hist.begin(); it2 != hist.end(); it2++) 
    {
        sorted_hist.push_back(make_pair(it2->first, it2->second));
    }

    // Sort the values based on the histogram count
    sort(sorted_hist.begin(), sorted_hist.end(), compare_int_desc);

    // Get the first 256 colors
    vector<pair<tuple<int, int, int>, int>> first_256;
    first_256 = vector<pair<tuple<int, int, int>, int>>(sorted_hist.begin(), sorted_hist.begin() + 256);
    return first_256;
}


bool compare_double_asc(const pair<tuple<int, int, int>, double> i, const pair<tuple<int, int, int>, double> j)
{
    return (i.second < j.second);
}

///////////////////////////////////////////////////////////////////////////////
//
//      Convert the image to an 8 bit image using populosity quantization.  
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Quant_Populosity()
{

    TargaImage::Quant_Uniform_32(data);

    vector<pair<tuple<int, int, int>, int>> first_256 = get_first_256_colors_by_popularity(height, width, data);

    // Go through every pixel and find the closest color from the colormap
    for(int i = 0; i < width * height * 4; i += 4)
    {
        int r = data[i + RED];
        int g = data[i + GREEN];
        int b = data[i + BLUE];

        // Compute the euclidean distance to the most popular colors
        vector<pair<tuple<int, int, int>, int>>::iterator it;

        double min_distance = 100000.0;
        int index = 0;
        int color_target_index = 0;

        for(it = first_256.begin(); it != first_256.end(); it++)
        {   
            double d_r = r - get<0>(it->first);
            double d_g = g - get<1>(it->first);
            double d_b = b - get<2>(it->first);
            double d = d_r*d_r + d_g*d_g + d_b*d_b;

            // If the distance was smaller than the previous min_distance, store the index
            if(d < min_distance)
            {
                min_distance = d;
                color_target_index = index;
            }
            index++;
        }

        // Insert the nearest color
        data[i + RED] = get<0>((first_256.begin() + color_target_index)->first);
        data[i + GREEN] = get<1>((first_256.begin() + color_target_index)->first);
        data[i + BLUE] = get<2>((first_256.begin() + color_target_index)->first);        
    }
    

    return true;
}// Quant_Populosity


///////////////////////////////////////////////////////////////////////////////
//
//      Dither the image using a threshold of 1/2.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Threshold()
{
    ClearToBlack();
    return false;
}// Dither_Threshold


///////////////////////////////////////////////////////////////////////////////
//
//      Dither image using random dithering.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Random()
{
    ClearToBlack();
    return false;
}// Dither_Random


///////////////////////////////////////////////////////////////////////////////
//
//      Perform Floyd-Steinberg dithering on the image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_FS()
{
    ClearToBlack();
    return false;
}// Dither_FS


///////////////////////////////////////////////////////////////////////////////
//
//      Dither the image while conserving the average brightness.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Bright()
{
    ClearToBlack();
    return false;
}// Dither_Bright


///////////////////////////////////////////////////////////////////////////////
//
//      Perform clustered differing of the image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Cluster()
{
    ClearToBlack();
    return false;
}// Dither_Cluster


///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using Floyd-Steinberg dithering over
//  a uniform quantization - the same quantization as in Quant_Uniform.
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Color()
{
    ClearToBlack();
    return false;
}// Dither_Color


///////////////////////////////////////////////////////////////////////////////
//
//      Composite the current image over the given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Over(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout <<  "Comp_Over: Images not the same size\n";
        return false;
    }

    ClearToBlack();
    return false;
}// Comp_Over


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "in" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_In(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_In: Images not the same size\n";
        return false;
    }

    ClearToBlack();
    return false;
}// Comp_In


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "out" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Out(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Out: Images not the same size\n";
        return false;
    }

    ClearToBlack();
    return false;
}// Comp_Out


///////////////////////////////////////////////////////////////////////////////
//
//      Composite current image "atop" given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Atop(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Atop: Images not the same size\n";
        return false;
    }

    ClearToBlack();
    return false;
}// Comp_Atop


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image with given image using exclusive or (XOR).  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Xor(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Xor: Images not the same size\n";
        return false;
    }

    ClearToBlack();
    return false;
}// Comp_Xor


///////////////////////////////////////////////////////////////////////////////
//
//      Calculate the difference bewteen this imag and the given one.  Image 
//  dimensions must be equal.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Difference(TargaImage* pImage)
{
    if (!pImage)
        return false;

    if (width != pImage->width || height != pImage->height)
    {
        cout << "Difference: Images not the same size\n";
        return false;
    }// if

    for (int i = 0 ; i < width * height * 4 ; i += 4)
    {
        unsigned char        rgb1[3];
        unsigned char        rgb2[3];

        RGBA_To_RGB(data + i, rgb1);
        RGBA_To_RGB(pImage->data + i, rgb2);

        data[i] = abs(rgb1[0] - rgb2[0]);
        data[i+1] = abs(rgb1[1] - rgb2[1]);
        data[i+2] = abs(rgb1[2] - rgb2[2]);
        data[i+3] = 255;
    }

    return true;
}// Difference


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 box filter on this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Box()
{
    double mask[5][5] = {{0.04, 0.04, 0.04, 0.04, 0.04}, 
                         {0.04, 0.04, 0.04, 0.04, 0.04},
                         {0.04, 0.04, 0.04, 0.04, 0.04},
                         {0.04, 0.04, 0.04, 0.04, 0.04},
                         {0.04, 0.04, 0.04, 0.04, 0.04}};
    filter_image_5x5(mask);
   
    return true;
}// Filter_Box


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Bartlett filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Bartlett()
{
    double mask[5][5] = {{0.0123, 0.0247, 0.0370, 0.0247, 0.0123}, 
                         {0.0247, 0.0494, 0.0741, 0.0494, 0.0247},
                         {0.0370, 0.0741, 0.1111, 0.0741, 0.0370},
                         {0.0247, 0.0494, 0.0741, 0.0494, 0.0247},
                         {0.0123, 0.0247, 0.0370, 0.0247, 0.0123}};
    filter_image_5x5(mask);
   
    return true;
}// Filter_Bartlett


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Gaussian()
{
    double mask[5][5] = {{0.0039, 0.0156, 0.0234, 0.0156, 0.0039},
						{0.0156, 0.0625, 0.0937, 0.0625, 0.0156},
						{0.0234, 0.0937, 0.1406, 0.0937, 0.0234},
						{0.0156, 0.0625, 0.0937, 0.0625, 0.0156},
						{0.0039, 0.0156, 0.0234, 0.0156, 0.0039}};
    filter_image_5x5(mask);
   
    return true;
}// Filter_Gaussian

///////////////////////////////////////////////////////////////////////////////
//
//      Perform NxN Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////

bool TargaImage::Filter_Gaussian_N( unsigned int N )
{
    // Dynamically allocate memory for the mask
    double **mask;
    mask = new double* [N];
    for(unsigned int i = 0; i < N; i++)
    {
        mask[i] = new double[N];
    }

    double denominator = pow(2.0, (2*N-2));
    
    double prev = 0;
    for(unsigned int i = 0; i < N; i++)
    {
        for(unsigned int j = 0; j < N; j++)
        {
            // Construct the pascal triangle with binomial
            double binomial = Binomial(N-2, j);
            double sum = binomial + prev;
            prev = binomial;
            if(i == 0)
            {
                mask[i][j] = sum * (1/denominator);
            } else
            {
                mask[i][j] = sum*mask[0][i];
            }
            
        }
    }

    filter_image_NxN(mask, N);

    // Free the mask memory
    for(unsigned int i = 0; i < N; i++)
    {
        delete[] mask[i];
    }
    delete[] mask;

    return true;
}// Filter_Gaussian_N


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 edge detect (high pass) filter on this image.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Edge()
{
    unsigned char *image_filtered_high_pass = new unsigned char[width * height * 4];
    memcpy(image_filtered_high_pass, data, sizeof(unsigned char)*width*height*4);

    // First we smooth the image "data"
    Filter_Gaussian();

    // Then subtract the smoothed image from the original, which equals to high pass filtering
    for(int i = 0; i < width*height; i++)
    {
        int r = image_filtered_high_pass[i*4 + RED] - data[i*4 + RED];
        int g = image_filtered_high_pass[i*4 + GREEN] - data[i*4 + GREEN];
        int b = image_filtered_high_pass[i*4 + BLUE] - data[i*4 + BLUE];
        // Set negative values to zero
        if(r < 0)
        {
           r = 0;
        } 

        if(g < 0)
        {
            g = 0;
        } 

        if(b < 0)
        {
            b = 0;
        } 

        image_filtered_high_pass[i*4 + RED] = r;
        image_filtered_high_pass[i*4 + GREEN] = g;
        image_filtered_high_pass[i*4 + BLUE] = b;
    }

    memcpy(data, image_filtered_high_pass, sizeof(unsigned char)*width*height*4);

    delete[] image_filtered_high_pass;
    return true;
}// Filter_Edge


///////////////////////////////////////////////////////////////////////////////
//
//      Perform a 5x5 enhancement filter to this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Enhance()
{
    unsigned char *edge_enchanced_image = new unsigned char[width * height * 4];
    memcpy(edge_enchanced_image, data, sizeof(unsigned char)*width*height*4);

    // Edge filter the image
    Filter_Edge();
    
    for(int i = 0; i < width*height; i++)
    {
        int r = edge_enchanced_image[i*4 + RED] + data[i*4 + RED];
        int g = edge_enchanced_image[i*4 + GREEN] + data[i*4 + GREEN];
        int b = edge_enchanced_image[i*4 + BLUE] + data[i*4 + BLUE];
        
        // Overflow
        if(r > 255)
        {
           r = 255;
        } 

        if(g > 255)
        {
            g = 255;
        } 

        if(b > 255)
        {
            b = 255;
        } 

        edge_enchanced_image[i*4 + RED] = r;
        edge_enchanced_image[i*4 + GREEN] = g;
        edge_enchanced_image[i*4 + BLUE] = b;
    }

    memcpy(data, edge_enchanced_image, sizeof(unsigned char)*width*height*4);

    delete[] edge_enchanced_image;
    return true;
}// Filter_Enhance


///////////////////////////////////////////////////////////////////////////////
//
//      Run simplified version of Hertzmann's painterly image filter.
//      You probably will want to use the Draw_Stroke funciton and the
//      Stroke class to help.
// Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::NPR_Paint()
{
    ClearToBlack();
    return false;
}



///////////////////////////////////////////////////////////////////////////////
//
//      Halve the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Half_Size()
{
    ClearToBlack();
    return false;
}// Half_Size


///////////////////////////////////////////////////////////////////////////////
//
//      Double the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Double_Size()
{
    ClearToBlack();
    return false;
}// Double_Size


///////////////////////////////////////////////////////////////////////////////
//
//      Scale the image dimensions by the given factor.  The given factor is 
//  assumed to be greater than one.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Resize(float scale)
{
    float w = 1 / scale;
    int x, y;

    for(int i = 0; i < width*height; i++)
    {
        x = i % width;
        y = i / width;
    }

    return false;
}// Resize


//////////////////////////////////////////////////////////////////////////////
//
//      Rotate the image clockwise by the given angle.  Do not resize the 
//  image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Rotate(float angleDegrees)
{
    ClearToBlack();
    return false;
}// Rotate


//////////////////////////////////////////////////////////////////////////////
//
//      Given a single RGBA pixel return, via the second argument, the RGB
//      equivalent composited with a black background.
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::RGBA_To_RGB(unsigned char *rgba, unsigned char *rgb)
{
    const unsigned char	BACKGROUND[3] = { 0, 0, 0 };

    unsigned char  alpha = rgba[3];

    if (alpha == 0)
    {
        rgb[0] = BACKGROUND[0];
        rgb[1] = BACKGROUND[1];
        rgb[2] = BACKGROUND[2];
    }
    else
    {
	    float	alpha_scale = (float)255 / (float)alpha;
	    int	val;
	    int	i;

	    for (i = 0 ; i < 3 ; i++)
	    {
	        val = (int)floor(rgba[i] * alpha_scale);
	        if (val < 0)
		    rgb[i] = 0;
	        else if (val > 255)
		    rgb[i] = 255;
	        else
		    rgb[i] = val;
	    }
    }
}// RGA_To_RGB


/* Uniform quantization of each color channel to 32 shades
*/
void TargaImage::Quant_Uniform_32(unsigned char *image)
{
    for(int h = 0; h < height; h++){
        int row_offset = h*width;
        
        for(int w = 0; w < width; w++){
            image[(row_offset + w)*4 + RED] = (unsigned char) floor(image[(row_offset + w)*4 + RED] / 8)*8;
            image[(row_offset + w)*4 + GREEN] = (unsigned char) floor(image[(row_offset + w)*4 + GREEN] / 8)*8;
            image[(row_offset + w)*4 + BLUE] = (unsigned char) floor(image[(row_offset + w)*4 + BLUE] / 8)*8;
        }
    }
}

void TargaImage::filter_image_5x5(double mask[][5])
{
    // For storing the average RGB values
    double avg_r, avg_g, avg_b;
    // For image navigation
    int x, y, pos_x, pos_y;
    int current_index;
    // Filter index
    int f_x, f_y;

    for(int i = 0; i < width*height; i++)
    {
        x = i % width;
        y = i / width;

        avg_r = 0;
        avg_g = 0;
        avg_b = 0;

        for(f_x = 0; f_x < 5; f_x++)
        {
            pos_x = x - 2 + f_x;

            if(pos_x < 0)
            {
                pos_x *= -1;
            } else if(pos_x > (width - 1))
            {
                pos_x = (width-1)*2 - pos_x;
            }

            for(f_y = 0; f_y < 5; f_y++)
            {
                pos_y = y - 2 + f_y;
                if(pos_y < 0)
                {
                    pos_y *= -1;
                } else if(pos_y > (height - 1))
                {
                    pos_y = (height-1)*2 - pos_y;
                }

                // Calculate the current position
                current_index = pos_x + pos_y*width;

                avg_r += data[current_index*4 + RED] * mask[f_x][f_y];
                avg_g += data[current_index*4 + GREEN] * mask[f_x][f_y];
                avg_b += data[current_index*4 + BLUE] * mask[f_x][f_y];
            }
        }
        data[i*4 + RED] = (int) avg_r;
        data[i*4 + GREEN] =  (int) avg_g;
        data[i*4 + BLUE] = (int) avg_b;
    }
}

void TargaImage::filter_image_NxN(double **mask, unsigned int mask_width)
{
    // For storing the average RGB values
    double avg_r, avg_g, avg_b;
    // For image navigation
    int x, y, pos_x, pos_y;
    int current_index;
    // Filter index
    unsigned int f_x, f_y;
    // Distance from filter center to edge
    int f_diff = (int) (mask_width / 2);
    
    for(int i = 0; i < width*height; i++)
    {
        x = i % width;
        y = i / width;

        avg_r = 0;
        avg_g = 0;
        avg_b = 0;

        for(f_x = 0; f_x < mask_width; f_x++)
        {
            pos_x = x - f_diff + f_x;

            if(pos_x < 0)
            {
                pos_x *= -1;
            } else if(pos_x > (width - 1))
            {
                pos_x = (width-1)*2 - pos_x;
            }

            for(f_y = 0; f_y < mask_width; f_y++)
            {
                pos_y = y - f_diff + f_y;
                if(pos_y < 0)
                {
                    pos_y *= -1;
                } else if(pos_y > (height - 1))
                {
                    pos_y = (height-1)*2 - pos_y;
                }

                // Calculate the current position
                current_index = pos_x + pos_y*width;

                avg_r += data[current_index*4 + RED] * mask[f_x][f_y];
                avg_g += data[current_index*4 + GREEN] * mask[f_x][f_y];
                avg_b += data[current_index*4 + BLUE] * mask[f_x][f_y];
            }
        }
        data[i*4 + RED] = (int) avg_r;
        data[i*4 + GREEN] =  (int) avg_g;
        data[i*4 + BLUE] = (int) avg_b;
    }

}

///////////////////////////////////////////////////////////////////////////////
//
//      Copy this into a new image, reversing the rows as it goes. A pointer
//  to the new image is returned.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage* TargaImage::Reverse_Rows(void)
{
    unsigned char   *dest = new unsigned char[width * height * 4];
    TargaImage	    *result;
    int 	        i, j;

    if (! data)
    	return NULL;

    for (i = 0 ; i < height ; i++)
    {
	    int in_offset = (height - i - 1) * width * 4;
	    int out_offset = i * width * 4;

	    for (j = 0 ; j < width ; j++)
        {
	        dest[out_offset + j * 4] = data[in_offset + j * 4];
	        dest[out_offset + j * 4 + 1] = data[in_offset + j * 4 + 1];
	        dest[out_offset + j * 4 + 2] = data[in_offset + j * 4 + 2];
	        dest[out_offset + j * 4 + 3] = data[in_offset + j * 4 + 3];
        }
    }

    result = new TargaImage(width, height, dest);
    delete[] dest;
    return result;
}// Reverse_Rows


///////////////////////////////////////////////////////////////////////////////
//
//      Clear the image to all black.
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::ClearToBlack()
{
    memset(data, 0, width * height * 4);
}// ClearToBlack


///////////////////////////////////////////////////////////////////////////////
//
//      Helper function for the painterly filter; paint a stroke at
// the given location
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::Paint_Stroke(const Stroke& s) {
   int radius_squared = (int)s.radius * (int)s.radius;
   for (int x_off = -((int)s.radius); x_off <= (int)s.radius; x_off++) {
      for (int y_off = -((int)s.radius); y_off <= (int)s.radius; y_off++) {
         int x_loc = (int)s.x + x_off;
         int y_loc = (int)s.y + y_off;
         // are we inside the circle, and inside the image?
         if ((x_loc >= 0 && x_loc < width && y_loc >= 0 && y_loc < height)) {
            int dist_squared = x_off * x_off + y_off * y_off;
            if (dist_squared <= radius_squared) {
               data[(y_loc * width + x_loc) * 4 + 0] = s.r;
               data[(y_loc * width + x_loc) * 4 + 1] = s.g;
               data[(y_loc * width + x_loc) * 4 + 2] = s.b;
               data[(y_loc * width + x_loc) * 4 + 3] = s.a;
            } else if (dist_squared == radius_squared + 1) {
               data[(y_loc * width + x_loc) * 4 + 0] = 
                  (data[(y_loc * width + x_loc) * 4 + 0] + s.r) / 2;
               data[(y_loc * width + x_loc) * 4 + 1] = 
                  (data[(y_loc * width + x_loc) * 4 + 1] + s.g) / 2;
               data[(y_loc * width + x_loc) * 4 + 2] = 
                  (data[(y_loc * width + x_loc) * 4 + 2] + s.b) / 2;
               data[(y_loc * width + x_loc) * 4 + 3] = 
                  (data[(y_loc * width + x_loc) * 4 + 3] + s.a) / 2;
            }
         }
      }
   }
}


///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke() {}

///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke(unsigned int iradius, unsigned int ix, unsigned int iy,
               unsigned char ir, unsigned char ig, unsigned char ib, unsigned char ia) :
   radius(iradius),x(ix),y(iy),r(ir),g(ig),b(ib),a(ia)
{
}

