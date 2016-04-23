//////////////////////////////////////////////////////////////////////
// BoundingBox.h: Interface
//
//////////////////////////////////////////////////////////////////////

#ifndef BOUNDING_BOX_H
#define BOUNDING_BOX_H

#include <TYPES.h>

#include <linearalgebra/Vector3.hpp>

#include <distancefield/distanceField.h>

#include <math.h>
#include <float.h>

using namespace std;

//////////////////////////////////////////////////////////////////////
// BoundingBox class
//
// Helpful class for storing object bounds
//////////////////////////////////////////////////////////////////////
class BoundingBox {
  public:
    BoundingBox( bool is3D = true )
      : _is3D( is3D )
    {
      _minBound = FLT_MAX;
      _maxBound = -FLT_MAX;
    }

    BoundingBox( const Vector3Array &vertices, bool is3D = true )
      : _is3D( is3D )
    {
      _minBound = FLT_MAX;
      _maxBound = -FLT_MAX;

      for ( size_t i = 0; i < vertices.size(); i++ )
      {
        if ( vertices.at( i )[0] > xmax() )
        {
          _maxBound[0] = vertices.at( i )[0];
        }
        if ( vertices.at( i )[1] > ymax() )
        {
          _maxBound[1] = vertices.at( i )[1];
        }
        
        if ( _is3D )
        {
          if ( vertices.at( i )[2] > zmax() )
          {
            _maxBound[2] = vertices.at( i )[2];
          }
        }

        if ( vertices.at( i )[0] < xmin() )
        {
          _minBound[0] = vertices.at( i )[0];
        }
        if ( vertices.at( i )[1] < ymin() )
        {
          _minBound[1] = vertices.at( i )[1];
        }

        if ( _is3D )
        {
          if ( vertices.at( i )[2] < zmin() )
          {
            _minBound[2] = vertices.at( i )[2];
          }
        }
      }
    }
   
    // make a squre bounding box centered at zero and have
    // cellSize*cellDivision length in each dimension
    BoundingBox(const REAL &cellSize, const int &cellDivisions,
                bool is3D=true)
        : _minBound(-cellSize*(REAL)cellDivisions/2.0,-cellSize*(REAL)cellDivisions/2.0,-cellSize*(REAL)cellDivisions/2.0), 
          _maxBound(+cellSize*(REAL)cellDivisions/2.0,+cellSize*(REAL)cellDivisions/2.0,+cellSize*(REAL)cellDivisions/2.0), 
          _is3D(is3D)
    {
    }

    BoundingBox( const Vector3d &minBound, const Vector3d &maxBound,
                 bool is3D = true )
      : _minBound( minBound ),
        _maxBound( maxBound ),
        _is3D( is3D )
    {
    }

    BoundingBox( const Vec3d &minBound, const Vec3d &maxBound,
                 bool is3D = true )
      : _is3D( is3D )
    {
      for ( int i = 0; i < 3; i++ )
      {
        _minBound[ i ] = minBound[ i ];
        _maxBound[ i ] = maxBound[ i ];
      }
    }

    BoundingBox( const FieldBoundingBox &bbox )
      : _is3D( true )
    {
      Vec3d bmin = bbox.bmin();
      Vec3d bmax = bbox.bmax();

      for ( int i = 0; i < 3; i++ )
      {
        _minBound[ i ] = bmin[ i ];
        _maxBound[ i ] = bmax[ i ];
      }
    }

    void             setMinBound( const Vector3d &minBound )
                     {
                       _minBound = minBound;
                     }
    void             setMaxBound( const Vector3d &maxBound )
                     {
                       _maxBound = maxBound;
                     }

    Vector3d        &minBound()
                     {
                         return _minBound;
                     }
    Vector3d        &maxBound()
                     {
                         return _maxBound;
                     }

    const Vector3d  &minBound() const
                     {
                         return _minBound;
                     }
    const Vector3d  &maxBound() const
                     {
                         return _maxBound;
                     }

    REAL             xmin() const
                     {
                         return _minBound[0];
                     }
    REAL             ymin() const
                     {
                        return _minBound[1];
                     }
    REAL             zmin() const
                     {
                        return _minBound[2];
                     }

    REAL             axismin( int i ) const
                     {
                        return _minBound[i];
                     }

    REAL             xmax() const
                     {
                        return _maxBound[0];
                     }
    REAL             ymax() const
                     {
                        return _maxBound[1];
                     }
    REAL             zmax() const
                     {
                        return _maxBound[2];
                     }

    REAL             axismax( int i ) const
                     {
                        return _maxBound[i];
                     }

    REAL             axislength( int i ) const
                     {
                        return axismax( i ) - axismin( i );
                     }
    REAL             maxlength() const
                     {
                        return max( axislength( 0 ), max( axislength( 1 ),
                                                         axislength( 2 ) ) );
                     }
    REAL             minlength() const 
                     {
                         return min( axislength(0), min(axislength(1), axislength(2))); 
                     }

    bool             isInside( const Vector3d &x ) const
                     {
                         if ( _is3D )
                         {
                             return ( x[0] >= xmin() &&
                                      x[1] >= ymin() &&
                                      x[2] >= zmin() &&
                                      x[0] <= xmax() &&
                                      x[1] <= ymax() &&
                                      x[2] <= zmax() );
                         }
                         else
                         {
                             return ( x[0] >= xmin() &&
                                      x[1] >= ymin() &&
                                      x[0] <= xmax() &&
                                      x[1] <= ymax() );
                         }
                     }

    // Returns the distance to the boundary of the box.
    // Negative if outside the box - positive if inside.
    REAL             boundaryDistance( const Vector3d &x ) const
                     {
                         REAL xdist = min( abs( x[0] - xmin() ),
                                           abs( x[0] - xmax() ) );
                         REAL ydist = min( abs( x[1] - ymin() ),
                                           abs( x[1] - ymax() ) );
                         REAL zdist = min( abs( x[2] - zmin() ),
                                           abs( x[2] - zmax() ) );
     
                         if ( isInside( x ) )
                             return min( xdist, min( ydist, zdist ) );
                         else
                             return -1.0 * min( xdist, min( ydist, zdist ) );
                     }

    Vector3d         center() const
                     {
                        return (REAL)0.5 * ( _minBound + _maxBound );
                     }

    // Whether or not two bounding boxes intersect
    bool             intersects( BoundingBox &b ) const
                     {
                        if ( _is3D )
                        {
                            return ( ( (xmin() <= b.xmax() && xmax() >= b.xmax())
                                    || (b.xmin() <= xmax() && b.xmax() >= xmax()) )
                                  && ( (ymin() <= b.ymax() && ymax() >= b.ymax())
                                    || (b.ymin() <= ymax() && b.ymax() >= ymax()) )
                                  && ( (zmin() <= b.zmax() && zmax() >= b.zmax())
                                    || (b.zmin() <= zmax() && b.zmax() >= zmax()) )
                                   );
                        }
                        else
                        {
                            return ( ( (xmin() <= b.xmax() && xmax() >= b.xmax())
                                     || (b.xmin() <= xmax() && b.xmax() >= xmax()) )
                                   && ( (ymin() <= b.ymax() && ymax() >= b.ymax())
                                     || (b.ymin() <= ymax() && b.ymax() >= ymax()) )
                                    );
                        }
                     }

    // Projects a point to the boundary of the box, if it lies outside.
    void             projectToBoundary( Vector3d &x ) const
                     {
                         x[0] = min( _maxBound[0], x[0] );
                         x[1] = min( _maxBound[1], x[1] );
     
                         if ( _is3D ) {
                             x[2] = min( _maxBound[2], x[2] );
                         }
     
                         x[0] = max( _minBound[0], x[0] );
                         x[1] = max( _minBound[1], x[1] );
     
                         if ( _is3D ) {
                             x[2] = max( _minBound[2], x[2] );
                         }
                     }

    // Intersects this bounding box with a vertex, making it
    // larger if necessary
    BoundingBox     &operator+=( const Vector3d &x )
                     {
                         if ( _is3D ) {
                             _minBound[0] = min( _minBound[0], x[0] );
                             _minBound[1] = min( _minBound[1], x[1] );
                             _minBound[2] = min( _minBound[2], x[2] );
         
                             _maxBound[0] = max( _maxBound[0], x[0] );
                             _maxBound[1] = max( _maxBound[1], x[1] );
                             _maxBound[2] = max( _maxBound[2], x[2] );
                         } else {
                             _minBound[0] = min( _minBound[0], x[0] );
                             _minBound[1] = min( _minBound[1], x[1] );
         
                             _maxBound[0] = max( _maxBound[0], x[0] );
                             _maxBound[1] = max( _maxBound[1], x[1] );
                         }

                         return *this;
                     }

    // Intersects this box with another one
    BoundingBox     &operator+=( BoundingBox &b )
                     {
                         if ( _is3D ) {
                             _minBound[0] = min( _minBound[0], b.xmin() );
                             _minBound[1] = min( _minBound[1], b.ymin() );
                             _minBound[2] = min( _minBound[2], b.zmin() );
         
                             _maxBound[0] = max( _maxBound[0], b.xmax() );
                             _maxBound[1] = max( _maxBound[1], b.ymax() );
                             _maxBound[2] = max( _maxBound[2], b.zmax() );
                         } else {
                             _minBound[0] = min( _minBound[0], b.xmin() );
                             _minBound[1] = min( _minBound[1], b.ymin() );

                             _maxBound[0] = max( _maxBound[0], b.xmax() );
                             _maxBound[1] = max( _maxBound[1], b.ymax() );
                         }

                         return *this;
                     }

    // Scales this bounding box, while keeping its center
    // in the same position
    BoundingBox     &operator*=( REAL alpha )
                     {
                         Vector3d minBound;
                         Vector3d maxBound;
     
                         minBound = alpha * ( _minBound - center() );
                         minBound += center();
     
                         maxBound = alpha * ( _maxBound - center() );
                         maxBound += center();
     
                         _minBound = minBound;
                         _maxBound = maxBound;
     
                         return *this;
                     }

  private:
    Vector3d           _minBound;
    Vector3d           _maxBound;

    bool               _is3D;
};

// Overloaded operators
inline BoundingBox operator*( BoundingBox &b, REAL alpha )
{
    Vector3d minBound = b.center() + alpha * ( b.minBound() - b.center() );
    Vector3d maxBound = b.center() + alpha * ( b.maxBound() - b.center() );

    return BoundingBox( minBound, maxBound );
}

#endif
