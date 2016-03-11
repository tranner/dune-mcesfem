/**************************************************************************
 
  The dune-fem module is a module of DUNE (see www.dune-project.org).
  It is based on the dune-grid interface library 
  extending the grid interface by a number of discretization algorithms
  for solving non-linear systems of partial differential equations.

  Copyright (C) 2003 - 2014 Robert Kloefkorn
  Copyright (C) 2003 - 2010 Mario Ohlberger 
  Copyright (C) 2004 - 2014 Andreas Dedner
  Copyright (C) 2005        Adrian Burri
  Copyright (C) 2005 - 2014 Mirko Kraenkel
  Copyright (C) 2006 - 2014 Christoph Gersbacher
  Copyright (C) 2006 - 2014 Martin Nolte
  Copyright (C) 2011 - 2014 Tobias Malkmus
  Copyright (C) 2012 - 2014 Stefan Girke
  Copyright (C) 2013 - 2014 Claus-Justus Heine
  Copyright (C) 2013 - 2014 Janick Gerstenberger
  Copyright (C) 2013        Sven Kaulman
  Copyright (C) 2013        Tom Ranner


  The dune-fem module is free software; you can redistribute it and/or 
  modify it under the terms of the GNU General Public License as 
  published by the Free Software Foundation; either version 2 of 
  the License, or (at your option) any later version.

  The dune-fem module is distributed in the hope that it will be useful, 
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along
  with this program; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 
**************************************************************************/
#ifndef HEAT_PROBLEMINTERFACE_HH
#define HEAT_PROBLEMINTERFACE_HH

#include <cassert>
#include <cmath>

#include <dune/fem/solver/timeprovider.hh>

#include "probleminterface.hh"

/** \brief problem interface class for time dependent problem descriptions, i.e. right hand side,
 *         boudnary data, and, if exsistent, an exact solution. A routine time() is
 *         provided. 
 */
template <class FunctionSpace>
class TemporalProblemInterface : public ProblemInterface<FunctionSpace>
{
public:
  typedef Dune::Fem::TimeProviderBase  TimeProviderType ;

  //! constructor taking time provider
  TemporalProblemInterface( const TimeProviderType &timeProvider ) 
    : timeProvider_(timeProvider)
  {
  }

  //! return current simulation time 
  double time() const
  {
    return timeProvider_.time();
  }

  //! return current time step size (\f$ \delta t \f$)
  double deltaT() const 
  {
    return timeProvider_.deltaT() ;
  }

  //! return reference to Problem's time provider 
  const TimeProviderType & timeProvider() const 
  {
    return timeProvider_;
  }

protected:
  const TimeProviderType &timeProvider_;
};
#endif // #ifndef HEAT_PROBLEMINTERFACE_HH

