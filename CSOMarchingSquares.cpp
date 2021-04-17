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
//!
/*!
// \file    CSOMarchingSquares.cpp
// \author  Alexander Koehn
// \date    03/2007
*/
//----------------------------------------------------------------------------------


#include "../CSOBase/CSO.h"
#include "CSOSmoothing.h"
#include "CSOMarchingSquares.h"
#include "CSOGeneratePathPoints.h"

ML_START_NAMESPACE

//////////////////////////////////////////////////////////////////////////

CSOMarchingSquares::CSOMarchingSquares()
{
  _image      = NULL;
  _imageSizeX = 0;
  _imageSizeY = 0;
  _startX     = 0;
  _startY     = 0;
  _isoValue   = 0;
  _function   = NULL;
  
  _voxelZPosition     = 0;  
  
  _bInterpolatePoints = false;
}

//////////////////////////////////////////////////////////////////////////

CSOMarchingSquares::~CSOMarchingSquares()
{
}

//////////////////////////////////////////////////////////////////////////

void CSOMarchingSquares::_createCell(int topLeftVoxel[2], CSOMarchingSquaresCell& cell)
{
  cell._isBorder = false;

  float values[4] = {0.f};

  const bool xBorderLeft   =  topLeftVoxel[0]    <   0;
  const bool xBorderRight  = (topLeftVoxel[0]+1) >= _imageSizeX;
  const bool yBorderTop    =  topLeftVoxel[1]    <   0;
  const bool yBorderBottom = (topLeftVoxel[1]+1) >= _imageSizeY;

  if (_image && (xBorderLeft || yBorderTop)) 
  {
    values[0] = _isoValue-1.f;
    cell._isBorder = true;
  } 
  else 
  {
    values[0] = _getValueAt(topLeftVoxel[0], topLeftVoxel[1]);
  }

  if (_image && (xBorderRight || yBorderTop)) 
  {
    values[1] = _isoValue-1.f;
    cell._isBorder = true;
  } 
  else 
  {
    values[1] = _getValueAt(topLeftVoxel[0]+1, topLeftVoxel[1]);
  }

  if (_image && (xBorderRight || yBorderBottom)) 
  {
    values[2] = _isoValue-1.f;
    cell._isBorder = true;
  } 
  else 
  {
    values[2] = _getValueAt(topLeftVoxel[0]+1, topLeftVoxel[1]+1);
  }

  if (_image && (xBorderLeft || yBorderBottom)) 
  {
    values[3] = _isoValue-1.f;
    cell._isBorder = true;
  } 
  else 
  {
    values[3] = _getValueAt(topLeftVoxel[0], topLeftVoxel[1]+1);
  }

  cell.set(topLeftVoxel, values, _isoValue);
}

//////////////////////////////////////////////////////////////////////////

void CSOMarchingSquares::_createCell(int topLeftVoxel[2], CSOMarchingSquaresCell& cell, float values[4], bool isBorder, char cellConfig)
{
  cell._isBorder = isBorder;
  cell.set(topLeftVoxel, values, _isoValue, cellConfig);
}

//////////////////////////////////////////////////////////////////////////

void CSOMarchingSquares::_createCell(int topLeftVoxel[2], std::set<unsigned int>& cellsVisited, std::queue<CSOMarchingSquaresCell>& cellFront)
{
  CSOMarchingSquaresCell cell;
  _createCell(topLeftVoxel, cell);
  const int key = _getKey(cell);
  if (cellsVisited.find(key) == cellsVisited.end())
  {
    cellFront.push(cell);
    cellsVisited.insert(key);
  }
}

//////////////////////////////////////////////////////////////////////////

int CSOMarchingSquares::_walkToCell(const CSOMarchingSquaresCell& fromCell, int fromDir, CSOMarchingSquaresCell& toCell)
{
  toCell.reset();

  // Get the direction in which to leave the fromCell, so we can determine the values for toCell.
  float values[4];
  const int   toDir     = fromCell.getToDirection(fromDir);
  int   topLeftVoxel[2] = { fromCell._topLeftVoxel[0], fromCell._topLeftVoxel[1] };

  // Use a different value for outside positions when evaluating a function
  const float borderIsoValueOffset = static_cast<float>((!_image && _function) ? -1 : 1);

  // Leaving the cell (fromCell) in direction toDir, enterDstDir will hold the
  // incoming direction for the next cell (toCell)
  int enterDstDir = CSOMarchingSquaresCell::NONE;

  switch(toDir)
  {
    case CSOMarchingSquaresCell::LEFT: // leaving fromDir to the left
      {
        enterDstDir = CSOMarchingSquaresCell::RIGHT;

        topLeftVoxel[0] -= 1;

        const bool xBorderLeft   =  topLeftVoxel[0] < _startX;
        const bool yBorderTop    =  topLeftVoxel[1] < _startY;
        const bool yBorderBottom = (topLeftVoxel[1]+1) >= _imageSizeY;

        if (xBorderLeft) 
        {
          values[0] = _isoValue-borderIsoValueOffset;
          values[3] = _isoValue-borderIsoValueOffset;
          toCell._isBorder = true;
        }

        if (yBorderTop) 
        {
          values[0] = _isoValue-borderIsoValueOffset;
          values[1] = _isoValue-borderIsoValueOffset;
          toCell._isBorder = true;
        } 
        else 
        {
          if (!xBorderLeft) { values[0] = _getValueAt(topLeftVoxel[0], topLeftVoxel[1]); } // _image[topLeftVoxel[0] + topLeftVoxel[1] * _imageSizeX]; }
          values[1] = fromCell._values[0];
        }

        if (yBorderBottom) 
        {
          values[3] = _isoValue-borderIsoValueOffset;
          values[2] = _isoValue-borderIsoValueOffset;
          toCell._isBorder = true;
        } 
        else 
        {
          if (!xBorderLeft) { values[3] = _getValueAt(topLeftVoxel[0], topLeftVoxel[1]+1); } // _image[topLeftVoxel[0] + (topLeftVoxel[1]+1) * _imageSizeX]; }
          values[2] = fromCell._values[3];
        }        
      } break;
    case CSOMarchingSquaresCell::TOP: // leaving fromDir to the top
      {
        enterDstDir = CSOMarchingSquaresCell::BOTTOM;
        topLeftVoxel[1] -= 1;

        const bool yBorderTop   =  topLeftVoxel[1] < _startY;
        const bool xBorderLeft  =  topLeftVoxel[0] < _startX;
        const bool xBorderRight = (topLeftVoxel[0]+1) >= _imageSizeX;

        if (yBorderTop) 
        {
          values[0] = _isoValue-borderIsoValueOffset;
          values[1] = _isoValue-borderIsoValueOffset;
          toCell._isBorder = true;
        }

        if (xBorderLeft) 
        {
          values[0] = _isoValue-borderIsoValueOffset;
          values[3] = _isoValue-borderIsoValueOffset;
          toCell._isBorder = true;
        } 
        else 
        {
          if (!yBorderTop) { values[0] = _getValueAt(topLeftVoxel[0], topLeftVoxel[1]); } // _image[topLeftVoxel[0] + topLeftVoxel[1] * _imageSizeX]; }
          values[3] = fromCell._values[0];
        }

        if (xBorderRight) 
        {
          values[1] = _isoValue-borderIsoValueOffset;
          values[2] = _isoValue-borderIsoValueOffset;
          toCell._isBorder = true;
        } 
        else 
        {
          if (!yBorderTop) { values[1] = _getValueAt(topLeftVoxel[0]+1, topLeftVoxel[1]); } // _image[topLeftVoxel[0] + 1 + topLeftVoxel[1] * _imageSizeX]; }
          values[2] = fromCell._values[1];
        }        
      } break;
    case CSOMarchingSquaresCell::RIGHT: // leaving fromDir to the right
      {
        enterDstDir = CSOMarchingSquaresCell::LEFT;
        topLeftVoxel[0] += 1;

        const bool xBorderRight  = (topLeftVoxel[0]+1) >= _imageSizeX;
        const bool yBorderTop    =  topLeftVoxel[1] < _startY;
        const bool yBorderBottom = (topLeftVoxel[1]+1) >= _imageSizeY;

        if (xBorderRight) 
        {
          values[1] = _isoValue-borderIsoValueOffset;
          values[2] = _isoValue-borderIsoValueOffset;
          toCell._isBorder = true;
        }

        if (yBorderTop) 
        {
          values[0] = _isoValue-borderIsoValueOffset;
          values[1] = _isoValue-borderIsoValueOffset;
          toCell._isBorder = true;
        } 
        else 
        {
          if (!xBorderRight) { values[1] = _getValueAt(topLeftVoxel[0]+1, topLeftVoxel[1]); } // _image[topLeftVoxel[0] + 1 + topLeftVoxel[1] * _imageSizeX]; }
          values[0] = fromCell._values[1];
        }

        if (yBorderBottom) 
        {
          values[3] = _isoValue-borderIsoValueOffset;
          values[2] = _isoValue-borderIsoValueOffset;
          toCell._isBorder = true;
        } 
        else 
        {
          if (!xBorderRight) { values[2] = _getValueAt(topLeftVoxel[0]+1, topLeftVoxel[1]+1); } // _image[topLeftVoxel[0] + 1 + (topLeftVoxel[1]+1) * _imageSizeX]; }
          values[3] = fromCell._values[2];
        }
      } break;
    case CSOMarchingSquaresCell::BOTTOM: // leaving fromDir to the bottom
      {
        enterDstDir = CSOMarchingSquaresCell::TOP;
        topLeftVoxel[1] += 1;

        const bool yBorderBottom = (topLeftVoxel[1]+1) >= _imageSizeY;
        const bool xBorderLeft   =  topLeftVoxel[0] < _startX;
        const bool xBorderRight  = (topLeftVoxel[0]+1) >= _imageSizeX;

        if (yBorderBottom)
        {
          values[2] = _isoValue-borderIsoValueOffset;
          values[3] = _isoValue-borderIsoValueOffset;
          toCell._isBorder = true;
        }

        if (xBorderLeft) 
        {
          values[0] = _isoValue-borderIsoValueOffset;
          values[3] = _isoValue-borderIsoValueOffset;
          toCell._isBorder = true;
        } 
        else 
        {
          if (!yBorderBottom) { values[3] = _getValueAt(topLeftVoxel[0], topLeftVoxel[1]+1); } // _image[topLeftVoxel[0] + (topLeftVoxel[1]+1) * _imageSizeX]; }
          values[0] = fromCell._values[3];
        }

        if (xBorderRight) 
        {
          values[1] = _isoValue-borderIsoValueOffset;
          values[2] = _isoValue-borderIsoValueOffset;
          toCell._isBorder = true;
        } 
        else 
        {
          if (!yBorderBottom) { values[2] = _getValueAt(topLeftVoxel[0]+1, topLeftVoxel[1]+1); } // _image[topLeftVoxel[0] + 1 + (topLeftVoxel[1]+1) * _imageSizeX]; }
          values[1] = fromCell._values[2];
        }
      } break;
    default:
      CSO_ERROR("Unknown direction in CSOMarchingSquares::_walkToCell. Should not happen!");
      break;
  }
  toCell.set(topLeftVoxel, values, _isoValue);

  return enterDstDir;
}

//////////////////////////////////////////////////////////////////////////

void CSOMarchingSquares::reset()
{
  _existingValues.clear();
}

//////////////////////////////////////////////////////////////////////////

void CSOMarchingSquares::findIsoLine(int startPosX, int startPosY, float isoValue, bool /*takeShortestPath*/, CSOMarchingSquaresCell::vecPoint2D& positions)
{
  if (!_image) 
  {
    CSO_ERROR("No image set in CSOMarchingSquares::findIsoLine.");
    return;
  }

  reset();

  _isoValue       = isoValue;
  int startPos[2] = { startPosX, startPosY };

  CSOMarchingSquaresCell startCell;
  if (_findNearestIsoCell(startPos, startCell)) 
  {
    int fromDir = _getPossibleEnterDirection(startCell);
    _trackIsoCell(startCell, fromDir, positions);
  }
}

//////////////////////////////////////////////////////////////////////////

void CSOMarchingSquares::findAllIsoLines(float isoValue, bool /*takeShortestPath*/, std::vector<CSOMarchingSquaresCell::vecPoint2D>& vecPositions)
{
  if (_image || _function)
  {
    reset();

    _isoValue = isoValue;
    int pos[2] = { 0 };

    // store visited cells in 2D char array, 1 == visited, 0 == not visited
    std::vector<char> visitedCells;
    size_t numCells = (_imageSizeX - _startX) * (_imageSizeY - _startY);
    visitedCells.resize(numCells);
    char* visitedPtr = &visitedCells[0];
    memset(visitedPtr, 0, numCells);

    CSOMarchingSquaresCell  cell;
    char* visitedPositionPtr = &visitedCells[0];

    if (_image) 
    {
      // Version for _image, which uses optimized image data access.
      float values[4];
      char cellConfig = 0;
      bool isBorder = false;
      bool isYBorder = false;
      const float borderValue = _isoValue - 1.f;
      float* imagePos0;
      float* imagePos1;
      for (pos[1] = _startY; pos[1] < _imageSizeY; ++pos[1])
      {
        imagePos0 = _image + _imageSizeX * pos[1];
        imagePos1 = imagePos0 + _imageSizeX;
        values[0] = *imagePos0++; //_getValueAt(_startX, pos[1]);
        if (pos[1] == _imageSizeY - 1) 
        {
          values[3] = borderValue;
          isYBorder = true;
        }
        else 
        {
          values[3] = *imagePos1++; // _getValueAt(_startX, pos[1] + 1);
        }

        isBorder = false;
        for (pos[0] = _startX; pos[0] < _imageSizeX; ++pos[0])
        {
          if (pos[0] != _imageSizeX - 1) 
          {
            values[1] = *imagePos0++; // _getValueAt(pos[0] + 1, pos[1]);
            if (!isYBorder) 
            {
              values[2] = *imagePos1++; // _getValueAt(pos[0] + 1, pos[1] + 1);
            }
            else 
            {
              isBorder = true;
              values[2] = borderValue;
            }
          }
          else 
          {
            isBorder = true;
            values[1] = borderValue;
            values[2] = borderValue;
          }
          // check if the value has been visited
          const char visited = *visitedPositionPtr++;
          if (visited == 0)
          {
            cellConfig = CSOMarchingSquaresCell::calculateCellConfig(_isoValue, values);
            if (CSOMarchingSquaresCell::isIsoCell(cellConfig)) 
            {
              _createCell(pos, cell, values, isBorder, cellConfig);
              const int fromDir = _getPossibleEnterDirection(cell);
              CSOMarchingSquaresCell::vecPoint2D positions;
              _trackIsoCell(cell, fromDir, positions, visitedPtr);
              vecPositions.push_back(positions);
            }
          }
          // move values to the left
          values[0] = values[1];
          values[3] = values[2];
        }
      }
    }
    else
    {
      // Old/unoptimized version for _function
      for (pos[1] = _startY; pos[1] < _imageSizeY; ++pos[1])
      {
        for (pos[0] = _startX; pos[0] < _imageSizeX; ++pos[0])
        {
          const char visited = *visitedPositionPtr++;
          if (visited == 0)
          {
            _createCell(pos, cell);
            if (cell.isIsoCell())
            {
              const int fromDir = _getPossibleEnterDirection(cell);
              CSOMarchingSquaresCell::vecPoint2D positions;
              _trackIsoCell(cell, fromDir, positions, visitedPtr);
              vecPositions.push_back(positions);
            }
          }
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////

bool CSOMarchingSquares::fillCSO(FillCSOParameters& parameters)
{
  CSO* cso = parameters.cso;

  if (!parameters.pImg || !cso || parameters.positions.empty())
  { 
    return false; 
  } 
  
  const unsigned int numPositions = static_cast<unsigned int>(parameters.positions.size());
  
  int smoothRange = parameters.smoothRange;
  if (parameters.useSmoothing)  
  {
    if (smoothRange < 1) { smoothRange = 1; }
    
    CSOSmoothing::smooth2DPositions(parameters.positions, numPositions, parameters.numSmoothPasses, parameters.smoothingFactor, smoothRange);
  }
  
  if (parameters.shouldReduceCSOToOnlyOneSeedPoint && parameters.smoothingMode == SMOOTHING_MODE_NONE)
  {
    // Optimized version for SMOOTHING_MODE_NONE and shouldReduceCSOToOnlyOneSeedPoint.
    // TODO: implement spline and other modes for this mode as well.
    cso->setInitialSeedAndPathPointsNoEvent(1);
    const Vector3 vc(parameters.positions[0][0], parameters.positions[0][1], parameters.posZ);
    Vector3 wc = parameters.pImg->mapVoxelToWorld(vc);
    cso->getSeedPointAt(0)->worldPosition = wc;

    // Instead of mapping each voxel position to world, we calculate the x/y world steps
    // and use the fact that the coordinates are in x,y on a slice only.
    CSOPathPoints* points = cso->getPathPointsAt(0);
    std::vector<Vector3> worldPositions;
    worldPositions.resize(numPositions);
    Vector3 offset = parameters.pImg->mapVoxelToWorld(Vector3(0,0,parameters.posZ));
    Vector3 xstep  = parameters.pImg->mapVoxelToWorld(Vector3(1, 0, parameters.posZ)) - offset;
    Vector3 ystep  = parameters.pImg->mapVoxelToWorld(Vector3(0, 1, parameters.posZ)) - offset;
    for (unsigned int i = 0; i < numPositions; ++i)
    {
      worldPositions[i] = offset + parameters.positions[i][0] * xstep + parameters.positions[i][1] * ystep;
    }
    worldPositions.push_back(worldPositions.at(0));
    points->setWorldPositions(worldPositions);
  }
  else
  {
    // Convert the scanned iso-line positions to seed points,
    // and generate and link path point lists.
    cso->setInitialSeedAndPathPointsNoEvent(numPositions);
    for (unsigned int i = 0; i < numPositions; ++i)
    {
      const Vector3 vc(parameters.positions[i][0], parameters.positions[i][1], parameters.posZ);
      Vector3 wc = parameters.pImg->mapVoxelToWorld(vc);
      cso->getSeedPointAt(i)->worldPosition = wc;
    }

    // Interpolate the path point lists
    switch (parameters.smoothingMode)
    {
    case SMOOTHING_MODE_NONE:
      for (unsigned int i = 0; i < cso->numPathPointLists(); ++i)
      {
        CSOGeneratePathPoints::fillPathPointsLinear(cso, cso->getPathPointsAt(i));
      }
      break;
    case SMOOTHING_MODE_SPLINE_APPROXIMATION:
      for (unsigned int i = 0; i < cso->numPathPointLists(); ++i)
      {
        CSOGeneratePathPoints::fillPathPointsSplineApproximation(cso, cso->getPathPointsAt(i));
      }
      break;
    case SMOOTHING_MODE_SPLINE_INTERPOLATION:
      for (unsigned int i = 0; i < cso->numPathPointLists(); ++i)
      {
        CSOGeneratePathPoints::fillPathPointsSplineInterpolation(cso, cso->getPathPointsAt(i));
      }
      break;
    default: break;
    }
    if (parameters.shouldReduceCSOToOnlyOneSeedPoint)
    {
      CSOGeneratePathPoints::reduceToOneSeedPoint(cso);
    }
  }
  
  return true;
}

//////////////////////////////////////////////////////////////////////////

void CSOMarchingSquares::_trackIsoCell(CSOMarchingSquaresCell startCell, int fromDir, CSOMarchingSquaresCell::vecPoint2D& positions, char* pVisitedCells)
{
  if (!startCell.isIsoCell()) 
  {
    return;
  }
  if (fromDir == 0) 
  {
    return;
  }

  int fromIdx    = 0;
  int toIdx      = 1;
  // The algorithm keeps track of the current cell and the next cell.
  CSOMarchingSquaresCell cell[2];
  cell[fromIdx] = startCell;

  const float xyOffset = (!_image && _function) ? 0.0f : 0.5f;

  do 
  {
    if (pVisitedCells) 
    {
      const int x = cell[fromIdx]._topLeftVoxel[0];
      const int y = cell[fromIdx]._topLeftVoxel[1];
      if ((x >= _startX) && (x < _imageSizeX) && (y >= _startY) && (y < _imageSizeY))
      {
        pVisitedCells[(x - _startX) + (y - _startY) * (_imageSizeX - _startX)] = 1;
      }
    }
    cell[fromIdx].addPoints(_bInterpolatePoints, fromDir, positions, xyOffset);

    fromDir = _walkToCell(cell[fromIdx], fromDir, cell[toIdx]);

    fromIdx = toIdx;
    toIdx   = !fromIdx;
  } while ((startCell != cell[fromIdx]) && (fromDir != 0));
}

//////////////////////////////////////////////////////////////////////////

bool CSOMarchingSquares::_findNearestIsoCell(int voxelPosX, int voxelPosY, CSOMarchingSquaresCell& cell)
{
  int voxelPos[2] = { voxelPosX, voxelPosY };
  return _findNearestIsoCell(voxelPos, cell);
}

//////////////////////////////////////////////////////////////////////////

bool CSOMarchingSquares::_findNearestIsoCell(int voxelPos[2], CSOMarchingSquaresCell& cell)
{
  bool bRet = false;

  CSOMarchingSquaresCell                  newCell;
  std::queue<CSOMarchingSquaresCell>      cellFront;
  std::set<unsigned int>  cellsVisited;

  _createCell(voxelPos, newCell);
  cellFront.push(newCell);
  int key = _getKey(newCell);
  cellsVisited.insert(key);

  while (!cellFront.empty()) 
  {
    CSOMarchingSquaresCell currCell = cellFront.front();
    cellFront.pop();

    // The cell is no iso cell -> create a new front out of that cell.
    if (!currCell.isIsoCell()) 
    {
      voxelPos[0] = currCell._topLeftVoxel[0]-1;
      voxelPos[1] = currCell._topLeftVoxel[1];
      _createCell(voxelPos, cellsVisited, cellFront);

      voxelPos[0] = currCell._topLeftVoxel[0] + 1;
      _createCell(voxelPos, cellsVisited, cellFront);
      
      voxelPos[0] = currCell._topLeftVoxel[0];
      voxelPos[1] = currCell._topLeftVoxel[1]-1;
      _createCell(voxelPos, cellsVisited, cellFront);

      voxelPos[1] = currCell._topLeftVoxel[1]+1;
      _createCell(voxelPos, cellsVisited, cellFront);      
    } 
    else 
    {
      cell = currCell;
      bRet = true;
      break;
    }
  }

  return bRet;
}

//////////////////////////////////////////////////////////////////////////

float CSOMarchingSquares::_getValueAt(int x, int y)
{
  float result = 0;
  
  if (_image)
  {   
	  if (_isInImage(x,y))
	  {
		  result = _image[x + y * _imageSizeX];
	  }
  } 
  else if (_function)
  { 
    const unsigned int key = _getKey(x, y);
    
    boost::unordered_map<int, double>::const_iterator existingPositionsIter = _existingValues.find(key);
    
    if (existingPositionsIter != _existingValues.end()) 
    {     
      result = static_cast<float>(existingPositionsIter->second);      
    } 
    else 
    {
      const Vector3 voxelPosition(x, y, _voxelZPosition);      
      result = _function->evaluateAtPos(voxelPosition);      
      _existingValues[key] = result;      
    }
  }
    
  return result;
}
 

//////////////////////////////////////////////////////////////////////////
bool CSOMarchingSquares::_isInImage(int x, int y) const
{
	return (x >= 0) && (y >= 0) && (x < _imageSizeX) && (y < _imageSizeY);
}

ML_END_NAMESPACE


