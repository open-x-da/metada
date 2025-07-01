#!/usr/bin/env python3
"""
Script to print WRF variables (ZNU, XLONG_U, XLAT_U, U) point by point
"""

import netCDF4 as nc
import numpy as np
import sys

def print_wrf_variables(filename):
    """Print WRF variables point by point"""
    
    # Open the NetCDF file
    try:
        dataset = nc.Dataset(filename, 'r')
    except Exception as e:
        print(f"Error opening file: {e}")
        return
    
    # Read variables
    try:
        znu = dataset.variables['ZNU'][:]      # (Time, bottom_top)
        xlong_u = dataset.variables['XLONG_U'][:]  # (Time, south_north, west_east_stag)
        xlat_u = dataset.variables['XLAT_U'][:]    # (Time, south_north, west_east_stag)
        u = dataset.variables['U'][:]          # (Time, bottom_top, south_north, west_east_stag)
        
        print(f"Variable dimensions:")
        print(f"ZNU: {znu.shape}")
        print(f"XLONG_U: {xlong_u.shape}")
        print(f"XLAT_U: {xlat_u.shape}")
        print(f"U: {u.shape}")
        print("-" * 80)
        
        # Get dimensions
        time_dim = u.shape[0]
        z_dim = u.shape[1] 
        y_dim = u.shape[2]
        x_dim = u.shape[3]
        
        # Print point by point (limiting output for readability)
        max_points = 20  # Limit to first 10 points for demonstration
        
        print(f"{'Time':<4} {'Z':<3} {'Y':<3} {'X':<3} {'ZNU':<8} {'XLONG_U':<10} {'XLAT_U':<10} {'U':<10}")
        print("-" * 80)
        
        point_count = 0
        for t in range(time_dim):
            for z in range(10, z_dim, 5):
                for y in range(20, y_dim, 20):
                    for x in range(30, x_dim, 20):
                        if point_count >= max_points:
                            print(f"... (showing only first {max_points} points)")
                            return
                            
                        znu_val = znu[t, z]
                        xlong_val = xlong_u[t, y, x]
                        xlat_val = xlat_u[t, y, x]
                        u_val = u[t, z, y, x]
                        
                        print(f"{t:<4} {z:<3} {y:<3} {x:<3} {znu_val:<8.4f} {xlong_val:<10.4f} {xlat_val:<10.4f} {u_val:<10.4f}")
                        point_count += 1
        
    except KeyError as e:
        print(f"Variable not found: {e}")
    except Exception as e:
        print(f"Error reading variables: {e}")
    finally:
        dataset.close()

def print_specific_point(filename, t=0, z=0, y=0, x=0):
    """Print specific point values"""
    
    try:
        dataset = nc.Dataset(filename, 'r')
        
        znu = dataset.variables['ZNU'][t, z]
        xlong_u = dataset.variables['XLONG_U'][t, y, x]
        xlat_u = dataset.variables['XLAT_U'][t, y, x]
        u = dataset.variables['U'][t, z, y, x]
        
        print(f"Point [t={t}, z={z}, y={y}, x={x}]:")
        print(f"  ZNU:     {znu:.6f}")
        print(f"  XLONG_U: {xlong_u:.6f}")
        print(f"  XLAT_U:  {xlat_u:.6f}")
        print(f"  U:       {u:.6f}")
        
        dataset.close()
        
    except Exception as e:
        print(f"Error: {e}")

def print_slice(filename, t=0, z=0, y_range=None, x_range=None):
    """Print a 2D slice of data"""
    
    try:
        dataset = nc.Dataset(filename, 'r')
        
        u_shape = dataset.variables['U'].shape
        print(f"Full U shape: {u_shape}")
        
        # Default ranges
        if y_range is None:
            y_range = (0, min(5, u_shape[2]))
        if x_range is None:
            x_range = (0, min(5, u_shape[3]))
            
        print(f"\nSlice at t={t}, z={z}, y={y_range[0]}:{y_range[1]}, x={x_range[0]}:{x_range[1]}")
        print('Y\\X'.ljust(3), end="")
        for x in range(x_range[0], x_range[1]):
            print(f"{x:>10}", end="")
        print()
        
        for y in range(y_range[0], y_range[1]):
            print(f"{y:<3}", end="")
            for x in range(x_range[0], x_range[1]):
                u_val = dataset.variables['U'][t, z, y, x]
                print(f"{u_val:>10.4f}", end="")
            print()
            
        dataset.close()
        
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python print_wrf_variables.py <wrfinput_file> [mode]")
        print("Modes:")
        print("  all     - Print first 10 points (default)")
        print("  point   - Print specific point (t=0,z=0,y=0,x=0)")
        print("  slice   - Print 2D slice (t=0,z=0,first 5x5 points)")
        sys.exit(1)
    
    filename = sys.argv[1]
    mode = sys.argv[2] if len(sys.argv) > 2 else "all"
    
    if mode == "all":
        print_wrf_variables(filename)
    elif mode == "point":
        print_specific_point(filename)
    elif mode == "slice":
        print_slice(filename)
    else:
        print(f"Unknown mode: {mode}")
