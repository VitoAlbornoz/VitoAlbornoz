// This function masks clouds in Landsat 8 imagery.
var maskClouds = function(image) {
  var scored = ee.Algorithms.Landsat.simpleCloudScore(image);
  return image.updateMask(scored.select(['cloud']).lt(20));
};

// This function masks clouds and adds quality bands to Landsat 8 images.
var addQualityBands = function(image) {
  return maskClouds(image)
    // NDVI
    .addBands(image.normalizedDifference(['B5', 'B4']))
    // time in days
    .addBands(image.metadata('system:time_start'));
};

// Load a 2014 Landsat 8 ImageCollection.
// Map the cloud masking and quality band function over the collection.
var collection = ee.ImageCollection('LANDSAT/LC08/C01/T1_TOA')
  .filterDate('2014-06-01', '2014-12-31')
  .map(addQualityBands);

// Create a cloud-free, most recent value composite.
var recentValueComposite = collection.qualityMosaic('system:time_start');

// Create a greenest pixel composite.
var greenestPixelComposite = collection.qualityMosaic('nd');

// Display the results.
Map.setCenter(-122.374, 37.8239, 12); // San Francisco Bay
var vizParams = {bands: ['B5', 'B4', 'B3'], min: 0, max: 0.4};
Map.addLayer(recentValueComposite, vizParams, 'recent value composite');
Map.addLayer(greenestPixelComposite, vizParams, 'greenest pixel composite');

// Compare to a cloudy image in the collection.
var cloudy = ee.Image('LANDSAT/LC08/C01/T1_TOA/LC08_044034_20140825');
Map.addLayer(cloudy, vizParams, 'cloudy');