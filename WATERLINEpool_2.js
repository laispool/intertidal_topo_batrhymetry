/////////////////////////////////////////////////////////////
// Defining the function to estimate topobathymetry        //
/////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////
///// Laís Pool  (lais.pool@gmail.com)                /////                        
///// Florianopolis, 05/06/2023                       /////              
///// Code Editor - Earth Engine                      /////                             
///////////////////////////////////////////////////////////
/* 
Da primeira etapa, foram extraídas imagens de NDWI cortadas pelos 
vetores. Na segunda etapa é preciso identificar o limite entre água 
e terra em cada uma delas

Coletar a coordenada de cada pixel e depois 
atribuir um valor de nível de maré correspondente (referência a um 
marégrafo). 
*/
////////////////////////////////////////////////////////////

var pathImage = 'projects/ee-index-images/assets/valeSand'; // change this to the path of the asset containing the NWI data
var pathTideg = 'projects/ee-tide-gauge/assets/'; // change this to the path of the asset containing the tide gauge data
var nameTideg = 'tideVale'; // name of the tide gauge observations file
var savePath = 'projects/ee-tide-gauge/assets/'; // change this to the path of the asset where you want to save the results
var tideGauge = ee.FeatureCollection(pathTideg + nameTideg); // loading tide gauge observations as a ee.FeatureCollection
print('First 10 elements of tide gauge:',tideGauge.limit(10));


// Define a function to combine the date and time columns into a new column
var combineDateTime = function(feature) {
  var date = ee.String(feature.get('Data'));
  var hour = ee.String(feature.get('Hora'));
  var combinedDateTime = ee.String(date).cat(' ').cat(hour);
  return feature.set('dateTime', combinedDateTime);
};
var updatedCollection = tideGauge.map(combineDateTime);

updatedCollection = updatedCollection.map(function(feature) {
  var dateTimeStr = ee.String(feature.get('dateTime')); 
  var dateTime = ee.Date.parse('dd/MM/yyyy HH:mm:ss', dateTimeStr); //dd/MM/yyyy HH:mm // yyyy-MM-dd HH:mm:ss  //dd/MM/yyyy HH:mm:ss
  return feature.set('dateTime', dateTime);
});
print('Updated FeatureCollection:', updatedCollection.limit(10));
tideGauge = updatedCollection;
//var ll = tideGauge.toList(tideGauge.size());
//print("data mais antiga", ll.get(-1));

/*Defining a function for mapping over the lists
var filterValues = function(elem) {
  return ee.Algorithms.If(elem[1].eq(1), elem[0], null);
};
*/      

// Defining the Otsu function
/*
declaration: function otsuThreshold(histogram):
This function was made to find the optimum threshold of a grey scaled image (single band)

    list of parametres:
    • histogram: an histogram found with:
        var histogram = image.reduceRegion({
          reducer: ee.Reducer.histogram(),
          geometry: geometry, 
          scale: 10
        });
        var imageHist = histogram.get('bandName');
      
      • image: the name of the dingle band image
      • bandName: the name of the band in your image
      
    global variables: (none)
    libraries needed: (none)
    return value: threshold value
*/
var otsuThreshold = function(histogram) {
  var counts = ee.Array(ee.Dictionary(imageHist).get('histogram')); 
  var means = ee.Array(ee.Dictionary(imageHist).get('bucketMeans')); 
  var size = means.length().get([0]);
  var total = counts.reduce(ee.Reducer.sum(), [0]).get([0]);
  var sum = means.multiply(counts).reduce(ee.Reducer.sum(), [0]).get([0]);
  var mean = sum.divide(total);

  var indices = ee.List.sequence(1, size);
  
  // Compute between sum of squares, where each mean partitions the data.
  var bss = indices.map(function(i) {
    var aCounts = counts.slice(0, 0, i);
    var aCount = aCounts.reduce(ee.Reducer.sum(), [0]).get([0]);
    var aMeans = means.slice(0, 0, i);
    var aMean = aMeans.multiply(aCounts)
        .reduce(ee.Reducer.sum(), [0]).get([0])
        .divide(aCount);
    var bCount = total.subtract(aCount);
    var bMean = sum.subtract(aCount.multiply(aMean)).divide(bCount);
    return aCount.multiply(aMean.subtract(mean).pow(2)).add(
          bCount.multiply(bMean.subtract(mean).pow(2)));
  });

  //print(ui.Chart.array.values(ee.Array(bss), 0, means));
  //print('the threshold is:    ', means.sort(bss).get([-1]));

  // Return the mean value corresponding to the maximum BSS.
  return means.sort(bss).get([-1]);
};

// Function to add datetime property to each image and correct the timezone
var addDatetimeProperty = function(image) {
  var id = ee.String(image.get('system:id'));  // Get the image ID
  var datetimeStr = id.slice(41, 56);  // Extract the datetime string from the ID 
  
  // Convert the datetime string to a datetime object
  var datetime = ee.Date.parse('yyyyMMdd\'T\'HHmmss', datetimeStr);
  var updatedDatetime = datetime.advance(-3, 'hour'); //for São Paulo Brazil: -3h timeZone
  // Return the image with the added datetime property
  return image.set('datetime', updatedDatetime);
};

// Combined function to calculate time difference and set properties
var calculateTimeDifferenceAndSetProperties = function(image, feature) {
  var imageDate = ee.Date(image.get('datetime'));
  var featureDate = ee.Date(feature.get('dateTime'));
  
  var diff = (featureDate.difference(imageDate, 'second')).abs();

  return feature.set({
    'timeDifference': ee.Number(diff),
    'dateImage': imageDate,
    'dateGauge': featureDate
  });
};

// Function to add the 'tide' value to each feature in a featureCollection
var addTideValue = function(feature) {
  var tideValue = closestFeatureImage.get('Dado');
  var dateImage = closestFeatureImage.get('dateImage');
  dateImage = ee.Date(dateImage).format('yyyy-MM-dd');
  var dateGauge = closestFeatureImage.get('dateGauge');
  dateGauge = ee.Date(dateGauge).format('yyyy-MM-dd');
  
  return feature.set({
    'Tide': tideValue,
    'dateImage': dateImage,
    'dateGauge': dateGauge
  });
};

// Define the main function
/*
declaration: function SDB_intertidal(pathImage, tideGauge, savePath):
This function finds the elevation of waterlines found using a threshold technique and tide gauge information

    list of parametres:
    • pathImage: nome of the asset path of the nwi images
    • tideGauge: name and path of the tide gauge file (on assets)
    • savePath: path where the results are going to be saved
      
    global variables: 
      • imageCollection: 
      • histogram: 
      • threshold:
    libraries needed: (none)
    return value: threshold value
*/

//var SDB_intertidal = function(pathImage, tg, savePath) {
  var assetList = ee.data.listImages(pathImage).images;
  var idsList = assetList.map(function(img) { return img.id}); 
  var imageCollection = ee.ImageCollection(idsList);
  var imageCollectionWithDatetime = imageCollection.map(addDatetimeProperty);
  print('Image Collection with Datetime:', imageCollectionWithDatetime);
  
  var numImages = imageCollectionWithDatetime.size();
  var mergedPointsWithTide = ee.FeatureCollection([]);
  var closestFeatureImage; 
  
  print('begining of problematic code');
  for (var i = 0; i < numImages.getInfo(); i++) {
    
    var nwi = ee.Image(imageCollectionWithDatetime.toList(numImages).get(i));
    Map.addLayer(nwi, {min:-1, max:1}, "NWI no mask",false);
     
    var hist_NDWI1 = ui.Chart.image.histogram({image:nwi.select('ndwi'), scale: 10, maxPixels: 10e9})
    .setOptions({title: 'Histograma  NWI sem corte '+i});
    print (hist_NDWI1);
    
    var histogram = nwi.reduceRegion({reducer: ee.Reducer.histogram(), tileScale: 16, maxPixels: 10e9});
    var imageHist = histogram.get('ndwi');
    var threshold = otsuThreshold(histogram.get('histogram'));

      var lessThreshold = nwi.lte(threshold);
      var filteredImage = nwi.updateMask(lessThreshold);
      Map.addLayer(filteredImage, {min:-1, max:1}, "NWI masked",false);
      var newImage = filteredImage;
      
      // Get the image's projection and scale
      var projection = newImage.projection();
      var scale = projection.nominalScale();
      
      var longitude = ee.Image.pixelLonLat().select('longitude');
      var latitude = ee.Image.pixelLonLat().select('latitude');
      var imageWithCoordinates = newImage.addBands(longitude.rename('longitude'))
                                   .addBands(latitude.rename('latitude')); // Add the longitude and latitude bands as properties to each pixel
      print("image With Coordinates", imageWithCoordinates);
      
      var hist_NDWI2 = ui.Chart.image.histogram({image:imageWithCoordinates.select('ndwi'), scale: 10, maxPixels: 10e9})
      .setOptions({title: 'Histograma  NWI lessThreshold '+i});
      print (hist_NDWI2);
      
      var percentage = 0.5; //15%
      print("thresholding da imagem " + i, threshold);

      var maskThreshold = newImage.expression(
        "(b >= (thrld - abs(thrld) * percentage)) && (b <= (thrld + abs(thrld) * percentage)) ? 1 : 0",  //
        {b: newImage.select('ndwi'), thrld: threshold, percentage: percentage});  
    
      var waterlineImage = imageWithCoordinates.updateMask(maskThreshold);
      print("waterlineImage", waterlineImage);
      Map.addLayer(waterlineImage.select('ndwi'), {}, "waterlineImage "+ i);
      
      var hist_NDWI3 = ui.Chart.image.histogram({image:waterlineImage.select('ndwi'), scale: 10, maxPixels: 10e9})
      .setOptions({title: 'Histograma  NWI com corte '+i});
      print (hist_NDWI3);
      
      //var maskNAN = waterlineImage.updateMask(waterlineImage.mask().not());
      var bounds = waterlineImage.geometry();
      Map.centerObject(bounds);
      
      var ndwiIntegral = waterlineImage.select('ndwi').int();
      var img2point = ndwiIntegral.reduceToVectors({
        geometry: bounds, 
        scale: 10,
        maxPixels: 10e9,
        tileScale: 16,
        geometryType: 'centroid', // Use centroid of each pixel as the feature geometry
        eightConnected: false // Consider only four-connected neighbors
      });
      print("img2point",img2point.limit(100));
      var featureList = img2point.toList(img2point.size());// Convert the feature collection back to a list
      print('featureList',featureList)
      var pointsCollection = ee.FeatureCollection(featureList);
      print('Points as feature:', pointsCollection.limit(100));
      
      // Calculate time difference and set properties for each feature
      var closestCollection = tideGauge.map(function(feature) {
        return calculateTimeDifferenceAndSetProperties(waterlineImage, feature);
      });

      
      var sortedCollection = closestCollection.sort('timeDifference');
      closestFeatureImage = sortedCollection.first(); // Update closestFeatureImage
  
      //print("data da imagem: ", closestFeatureImage.get('dateImage'))
      //print("data da maré",closestFeatureImage.get('dateTime'));
      
      var tideValue = closestFeatureImage.get('Dado');
      //print("dado mais perto", tideValue)

      var pointsWithTide = pointsCollection.map(addTideValue);
      //print("pointsWithTide", pointsWithTide);
      
      // Add the pointsWithTide FeatureCollection to the list
      mergedPointsWithTide = mergedPointsWithTide.merge(pointsWithTide);
  }

  
//print("Merged Points with Tide:", mergedPointsWithTide.limit(10));

// Export the merged FeatureCollection
Export.table.toAsset({
  collection: mergedPointsWithTide,
  description: 'pointsWithTide_exportAsset',
  assetId: 'projects/ee-tide-gauge/assets/pointsWithTide_export'
});

Export.table.toDrive({
  collection: mergedPointsWithTide,
  description: 'pointsWithTide_exportCSV',
  fileFormat: 'CSV',
  folder:'Intertidal_Zones'
});

Export.table.toDrive({
  collection: mergedPointsWithTide,
  description: 'pointsWithTide_export',
  fileFormat: 'SHP',
  folder:'Intertidal_Zones'
});

