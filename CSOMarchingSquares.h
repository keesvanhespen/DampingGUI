// -----------------------------------------------------------------------
// 
// Copyright (c) 2001-2018, MeVis Medical Solutions AG, Bremen, Germany
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of MeVis Medical Solutions AG nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY MEVIS MEDICAL SOLUTIONS AG ''AS IS'' AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL MEVIS MEDICAL SOLUTIONS AG BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//----------------------------------------------------------------------------------
//! This class implements the Marching Squares algorithm to find isolines on 2D image slices.
/*!
// \file    CSOMarchingSquares.h
// \author  Alexander Koehn
// \date    03/2007
*/
//----------------------------------------------------------------------------------

#pragma once


#include "MLCSOIncludes.h"
#include "CSOMarchingSquaresCell.h"
#include "CSOFunction.h"

#include <ThirdPartyWarningsDisable.h>
#include <boost/unordered_map.hpp>
#include <ThirdPartyWarningsRestore.h>

#include <queue>

ML_START_NAMESPACE


//////////////////////////////////////////////////////////////////////////

//! This class implements the Marching Squares algorithm to find isolines on 2D image slices.
//! One can either find an iso line starting from a start point or find all contours on that slice.
class MLCSO_EXPORT CSOMarchingSquares
{
public:

  //! Constructor.
  CSOMarchingSquares();
  //! Destructor.
  ~CSOMarchingSquares();

  
public:

  //! Resets the state of the object.
  void reset();
  //! Specifies whether the algorithm should interpolate the contour bi-linearly.
  void setInterpolation(bool use);
  //! Sets the 2D image to track/find contours on.
  void setImage(float* image, int imgSizeX, int imgSizeY);
  //! Sets a 3D interpolation function, an image size and voxel size and a fixed z-offset.
  //! NOTE that this approach supports an axial scanning of the function only at the moment;
  //!      if this class should scan in sagittal or coronal direction as well, some extra flag needs to be provided.
  void setFunction(CSOFunction* implicitFunction, const Matrix4& voxelToWorldMatrix, int startX, int startY, int imgSizeX, int imgSizeY, int voxelZPosition);
  
  //! Looks for an isoline on the \c image. It starts searching for a starting cell using the startPos parameters and
  //! then tracks this cell.
  //! All contour points are pushed back to the positions vector.
  void findIsoLine(int startPosX, int startPosY, float isoValue, bool takeShortestPath, CSOMarchingSquaresCell::vecPoint2D& positions);
  //! Finds all contours on the \c image with isovalue \c isovalue. Returns them in a vector of 2D position vectors.
  void findAllIsoLines(float isoValue, bool takeShortestPath, std::vector<CSOMarchingSquaresCell::vecPoint2D>& vecPositions);

  struct FillCSOParameters
  {
    CSO* cso;
    CSOMarchingSquaresCell::vecPoint2D positions;
    float posZ;
    PagedImage* pImg;
    bool useSmoothing;
    float smoothingFactor;
    int numSmoothPasses;
    int smoothRange;
    CSOSmoothingModes smoothingMode;
    bool shouldReduceCSOToOnlyOneSeedPoint;
  };

  //! Fills the given CSO with seed points according to the given list of points.
  bool fillCSO(FillCSOParameters& parameters);
  
protected:

  //! Searches for a starting position such that startX/startY make up the top-left voxel position of a cell which is intersected by the isoline.
  void _findStartPosition(int& startX, int& startY);
  //! Searches for the first cell that is an iso cell starting from [\c voxelPosX, \c voxelPosY]. Returns true if a cell was found, i.e. \c cell is valid.
  bool _findNearestIsoCell(int voxelPosX, int voxelPosY, CSOMarchingSquaresCell& cell);
  //! Searches for the first cell that is an iso cell starting from [\c voxelPos[0], \c voxelPos[1]]. Returns true if a cell was found, i.e. \c cell is valid.
  bool _findNearestIsoCell(int voxelPos[2], CSOMarchingSquaresCell& cell);
  //! Tracks the isoline starting from the \c startCell. If this is not an iso-cell (i.e. the isoline intersects this cell) it just returns.
  //! One can optionally pass a map where to save all visited cells.
  void _trackIsoCell(CSOMarchingSquaresCell startCell, int fromDir, CSOMarchingSquaresCell::vecPoint2D& positions, char* pVisitedCells = NULL);
  //! Create a new cell at position \c topLeftVoxel[2]
  void _createCell(int topLeftVoxel[2], CSOMarchingSquaresCell& cell);
  //! Create a new cell at position \c topLeftVoxel[2], keeps track whether this cell has already been created.
  void _createCell(int topLeftVoxel[2], std::set<unsigned int>& cellsVisited, std::queue<CSOMarchingSquaresCell>& cellFront);
  //! Create a new cell at position \c topLeftVoxel[2]
  void _createCell(int topLeftVoxel[2], CSOMarchingSquaresCell& cell, float values[4], bool isBorder, char cellConfig);
  //! Create a cell \c toCell from the cell \c fromCell coming from \c fromDir.
  //! So this method first looks for the direction leaving the fromCell coming from the direction \c fromDir.
  //! @return The direction from which one enters the \c toCell
  int  _walkToCell(const CSOMarchingSquaresCell& fromCell, int fromDir, CSOMarchingSquaresCell& toCell);
  //! @return A start direction from which one can enter the created cell.
  int _getPossibleEnterDirection(const CSOMarchingSquaresCell& cell) const;
  //! Returns the key for a cell. E.g. for hashtables.
  unsigned int  _getKey(const CSOMarchingSquaresCell& cell) const;
  //! Returns the key for a position for hashtables.
  unsigned int _getKey(int x, int y) const;
  
  //! Returns the value of the image or function at the given position.
  float _getValueAt(int x, int y);

  bool _isInImage(int x, int y) const;
  
protected:

  //! Pointer to input image.
  float* _image;
  //! The x extent of the input image.
  int _imageSizeX;
  //! The y extent of the input image.
  int _imageSizeY;
  
  //! Starting voxel x (function).
  int _startX;
  //! Starting voxel y (function).
  int _startY;
  
  //! The iso value to find the iso line for.
  float _isoValue;
  //! An implicit function.
  CSOFunction* _function;
  //! The z-position of the voxels.
  int _voxelZPosition;
  //! The voxelToWorld matrix.
  Matrix4 _voxelToWorldMatrix;
  
  //! Should the contour be interpolated bi-linearly?
  bool _bInterpolatePoints;

  //! A map holding all computed image/function values.
  boost::unordered_map <int, double> _existingValues;
};

//////////////////////////////////////////////////////////////////////////

inline void CSOMarchingSquares::setImage(float* image, int imgSizeX, int imgSizeY)
{
  _image = image;
  _imageSizeX = imgSizeX;
  _imageSizeY = imgSizeY;
  _startX = 0;
  _startY = 0;
  
  _function = NULL;
  _voxelZPosition = 0; 
}

//////////////////////////////////////////////////////////////////////////

inline void CSOMarchingSquares::setFunction(CSOFunction* implicitFunction,
                                            const Matrix4& voxelToWorldMatrix,
                                            int startX, int startY,
                                            int imgSizeX, int imgSizeY, 
                                            int voxelZPosition)
{
  _image      = NULL;
  
  _imageSizeX = imgSizeX;
  _imageSizeY = imgSizeY;
  
  _startX = startX;
  _startY = startY;
  
  _function           = implicitFunction;
  _voxelToWorldMatrix = voxelToWorldMatrix;
  
  _voxelZPosition     = voxelZPosition;    
}

//////////////////////////////////////////////////////////////////////////

inline void CSOMarchingSquares::setInterpolation(bool useInterpolation)
{
  _bInterpolatePoints = useInterpolation;
}

//////////////////////////////////////////////////////////////////////////

inline unsigned int CSOMarchingSquares::_getKey(const CSOMarchingSquaresCell& cell) const
{
  return _getKey(cell._topLeftVoxel[0], cell._topLeftVoxel[1]);
}

//////////////////////////////////////////////////////////////////////////

inline unsigned int CSOMarchingSquares::_getKey(int x, int y) const
{
  return 1 + x-_startX + ((1 + y-_startY) * (_imageSizeX+_startX+1));
}

//////////////////////////////////////////////////////////////////////////

inline int CSOMarchingSquares::_getPossibleEnterDirection(const CSOMarchingSquaresCell& cell) const
{
  int fromDir = cell.getToDirection(1);
  if (0 == fromDir) {
    fromDir = cell.getToDirection(2);
    if (0 == fromDir) {
      fromDir = cell.getToDirection(4);
      if (0 == fromDir) {
        fromDir = cell.getToDirection(8);
      }
    }
  }
  return fromDir;
}

//////////////////////////////////////////////////////////////////////////

ML_END_NAMESPACE
