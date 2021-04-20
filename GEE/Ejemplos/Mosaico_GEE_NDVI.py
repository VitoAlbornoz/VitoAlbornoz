# %%
"""
<table class="ee-notebook-buttons" align="left">
    <td><a target="_blank"  href="https://github.com/giswqs/geemap/tree/master/examples/template/template.ipynb"><img width=32px src="https://www.tensorflow.org/images/GitHub-Mark-32px.png" /> View source on GitHub</a></td>
    <td><a target="_blank"  href="https://nbviewer.jupyter.org/github/giswqs/geemap/blob/master/examples/template/template.ipynb"><img width=26px src="https://upload.wikimedia.org/wikipedia/commons/thumb/3/38/Jupyter_logo.svg/883px-Jupyter_logo.svg.png" />Notebook Viewer</a></td>
    <td><a target="_blank"  href="https://colab.research.google.com/github/giswqs/geemap/blob/master/examples/template/template.ipynb"><img src="https://www.tensorflow.org/images/colab_logo_32px.png" /> Run in Google Colab</a></td>
</table>
"""

# %%
"""
## Install Earth Engine API and geemap
Install the [Earth Engine Python API](https://developers.google.com/earth-engine/python_install) and [geemap](https://geemap.org). The **geemap** Python package is built upon the [ipyleaflet](https://github.com/jupyter-widgets/ipyleaflet) and [folium](https://github.com/python-visualization/folium) packages and implements several methods for interacting with Earth Engine data layers, such as `Map.addLayer()`, `Map.setCenter()`, and `Map.centerObject()`.
The following script checks if the geemap package has been installed. If not, it will install geemap, which automatically installs its [dependencies](https://github.com/giswqs/geemap#dependencies), including earthengine-api, folium, and ipyleaflet.
"""

# %%
# Installs geemap package
import subprocess

try:
    import geemap
except ImportError:
    print("Installing geemap ...")
    subprocess.check_call(["python", "-m", "pip", "install", "geemap"])

# %%
import ee
import geemap

# %%
"""
## Create an interactive map 
The default basemap is `Google Maps`. [Additional basemaps](https://github.com/giswqs/geemap/blob/master/geemap/basemaps.py) can be added using the `Map.add_basemap()` function. 
"""

# %%
Map = geemap.Map(center=[40, -100], zoom=4)
Map

# %%
"""
## Add Earth Engine Python script 
"""

# %%
# Add Earth Engine dataset
# This function masks clouds in Landsat 8 imagery.
def maskClouds(image):
  scored = ee.Algorithms.Landsat.simpleCloudScore(image)
  return image.updateMask(scored.select(['cloud']).lt(20))


# This function masks clouds and adds quality bands to Landsat 8 images.
def addQualityBands(image):
  return maskClouds(image) \
    .addBands(image.normalizedDifference(['B5', 'B4'])) \
    .addBands(image.metadata('system:time_start'))


# Load a 2014 Landsat 8 ImageCollection.
# Map the cloud masking and quality band function over the collection.
collection = ee.ImageCollection('LANDSAT/LC08/C01/T1_TOA') \
  .filterDate('2014-06-01', '2014-12-31') \
  .map(addQualityBands)

# Create a cloud-free, most recent value composite.
recentValueComposite = collection.qualityMosaic('system:time_start')

# Create a greenest pixel composite.
greenestPixelComposite = collection.qualityMosaic('nd')

# Display the results.
Map.setCenter(-122.374, 37.8239, 12); # San Francisco Bay
vizParams = {'bands': ['B5',  'B4',  'B3'], 'min': 0, 'max': 0.4}
Map.addLayer(recentValueComposite, vizParams, 'recent value composite')
Map.addLayer(greenestPixelComposite, vizParams, 'greenest pixel composite')

# Compare to a cloudy image in the collection.
cloudy = ee.Image('LANDSAT/LC08/C01/T1_TOA/LC08_044034_20140825')
Map.addLayer(cloudy, vizParams, 'cloudy')


# %%
"""
## Display Earth Engine data layers 
"""

# %%
Map.addLayerControl()  # This line is not needed for ipyleaflet-based Map.
Map
