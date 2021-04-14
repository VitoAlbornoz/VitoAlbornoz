// Example 1 – IMAGE VISUALIZATION AND METADATA
// More info at: https://developers.google.com/earth-engine/image_visualization

// Load a LANDSAT 8 image from San Francisco Bay Area
var image = ee.Image('LANDSAT/LC8_L1T_TOA/LC80440342014077LGN00');

// Print image metadata
print(image);

// Define TRUE color visualization parameters (red -> R, green -> G, blue -> B)
var visualParTrue = {
    bands: ['B4', 'B3', 'B2'],
    min: 0,
    max: 0.2
};

// Define FALSE color visualization parameters #1 (red -> NIR, green -> R, blue -> G)
var visualParFalse1 = {
    bands: ['B5', 'B4', 'B3'],
    min: 0,
    max: 0.5
};

// Define FALSE color visualization parameters #2 (red -> G, green -> NIR, blue -> UB)
var visualParFalse2 = {
    bands: ['B3', 'B5', 'B1'],
    min: 0,
    max: 0.5
};

// Center map on the image and display it in true and false colors
Map.centerObject(image, 8);
Map.addLayer(image, visualParTrue, 'true color composite', true);
Map.addLayer(image, visualParFalse1, 'false color composite #1', false);
Map.addLayer(image, visualParFalse2, 'false color composite #2', false);

/*
Info: if its hard to visualize the images, you can manualy set the gamma correction factors using the Map interface
1) go to -> layer manager (Layers)
2) hover over one of the layers
3) left click on the cogwheel symbol
4) adjust the Gamma bar
5) left click on Apply button
*/

