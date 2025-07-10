# Simple Backend Tutorial Data

## Data Source Attribution

The tutorial datasets in this directory are **copied and adapted** from the [PDAF (Parallel Data Assimilation Framework)](http://pdaf.awi.de/) tutorial input datasets.

### Original Source
- **Project**: PDAF - Parallel Data Assimilation Framework
- **Website**: http://pdaf.awi.de/
- **License**: LGPL (Lesser General Public License)
- **Authors**: Lars Nerger, Martin Drews, and PDAF development team

### Adaptations Made
The original PDAF tutorial data has been adapted for use with the Metada framework:
- File formats may have been converted or modified for compatibility
- Data structure may have been adjusted to match Metada's simple backend interface
- Configuration files (YAML/JSON) have been created specifically for Metada

### Files Description

#### Ensemble Data
- `ens_1.txt` to `ens_9.txt`: Ensemble member state files
- Each file contains gridded state data for data assimilation experiments

#### Observation Data
- `obs.txt`: Primary observation dataset
- `obsB.txt`: Alternative observation dataset  
- `obsC.txt`: Compact observation dataset

#### Configuration Files
- `etkf.yaml` / `etkf.json`: Configuration for Ensemble Transform Kalman Filter
- `letkf.yaml` / `letkf.json`: Configuration for Local Ensemble Transform Kalman Filter

### Usage
These files serve as tutorial examples for testing and demonstrating the Metada framework's data assimilation capabilities with the simple backend.

### Acknowledgments
We gratefully acknowledge the PDAF development team for providing high-quality tutorial datasets that serve as the foundation for these examples. The original PDAF tutorial materials have been invaluable for developing and testing data assimilation algorithms.

### License
The adapted data retains compatibility with the original PDAF licensing terms. Please refer to the PDAF project for specific license details regarding the original datasets. 