#include <set>
#include "image.hpp"


double image_compare(image &trueimg, image &infered) {
    assert(trueimg.rows() == infered.rows());
    assert(trueimg.cols() == infered.cols());
    // get the set of colors in the trueimg
    std::set<int> colors;
    for (size_t i = 0; i < trueimg.rows(); ++i) {
      for (size_t j = 0; j < trueimg.cols(); ++j) {
        colors.insert(size_t(trueimg.pixel(i,j)));
      }
    }
    
    // fill a rounding color map
    int colormap[256];
    int previval = -256;
    std::set<int>::iterator curi = colors.begin();
    std::set<int>::iterator nexti = curi;
    nexti++;
    int nextival = (nexti != colors.end())?*nexti:512;
    while (curi != colors.end()) {
      int low = (previval + (*curi)) / 2; if (low < 0) low = 0;
      int high = (nextival + (*curi)) / 2; if (high > 256) high = 256;
      
      for (int i = low; i < high; ++i) {
          colormap[i] = (*curi);
      }
      previval = (*curi);
      curi++;
      nexti++;
      nextival = (nexti != colors.end())?*nexti:512;
    }
    for (size_t i = 0; i < 256; ++i) std::cout << colormap[i] << " ";
    std::cout << std::endl;
    // round the infered image
    for (size_t i = 0; i < infered.rows(); ++i) {
      for (size_t j = 0; j < infered.cols(); ++j) {
        if (infered.pixel(i,j) >= 255) infered.pixel(i,j) = 255;
        if (infered.pixel(i,j) < 0) infered.pixel(i,j) = 0;
        infered.pixel(i,j) = colormap[(size_t)(infered.pixel(i,j) + 0.5)];
      }
    }
    
    // compute difference
    double rmse = 0;
    for (size_t i = 0; i < infered.rows(); ++i) {
      for (size_t j = 0; j < infered.cols(); ++j) {
        rmse += (infered.pixel(i,j) - trueimg.pixel(i,j))*
                      (infered.pixel(i,j) - trueimg.pixel(i,j));
      }
    }
    rmse /= (infered.rows() * infered.cols());
    rmse = sqrt(rmse);
    return rmse;
}