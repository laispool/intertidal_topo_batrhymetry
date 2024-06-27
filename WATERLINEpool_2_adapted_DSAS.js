/////////////////////////////////////////////////////////////
// Defining the function to estimate topobathymetry        //
/////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////
///// Laís Pool  (lais.pool@gmail.com)                /////                        
///// Florianopolis, 05/06/2023                       /////              
///// Code Editor - Earth Engine                      /////                             
///////////////////////////////////////////////////////////
/* 
Algoritmo para inferir valores de topobatimetria na zona de intermaré.
Retirado de Costa etal 2021.

Da primeira etapa, foram extraídas imagens de NDWI cortadas pelos 
vetores. Na segunda etapa é preciso identificar o limite entre água 
e terra em cada uma delas

*/
////////////////////////////////////////////////////////////

var pathImage = 'projects/ee-index-images/assets/delta'; // change this to the path of the asset containing the NWI data
var pathTideg = 'projects/ee-tide-gauge/assets/'; // change this to the path of the asset containing the tide gauge data
var nameTideg = 'tide_Itapoa_delta'; // name of the tide gauge observations file
var savePath = 'projects/ee-tide-gauge/assets/'; // change this to the path of the asset where you want to save the results
var tideGauge = ee.FeatureCollection(pathTideg + nameTideg); // loading tide gauge observations as a ee.FeatureCollection
print('First 10 elements of tide gauge:',tideGauge.limit(10));

// algumas definições
var scale = 10;
var geometry = ee.Geometry.Polygon(
[[[-48.6064338684082,-26.22352298843246],
                [-48.55304718017578,-26.22352298843246],
                [-48.55304718017578,-26.185018250078308],
                [-48.6064338684082,-26.185018250078308],
                [-48.6064338684082,-26.22352298843246]]]);
Map.centerObject(geometry, 14);
Map.addLayer(geometry, {'color': 'green' }, "ROI", false);



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
  var dateTime = ee.Date.parse('dd/MM/yyyy HH:mm', dateTimeStr); //yyyy-MM-dd HH:mm:ss 
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
var thresholdingAlgorithm = function(histogram) { //otsuThreshold
  var counts = ee.Array(ee.Dictionary(histNDWI).get('histogram')); 
  var means = ee.Array(ee.Dictionary(histNDWI).get('bucketMeans')); 
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
  var datetimeStr = id.slice(38, 53);  // Extract the datetime string from the ID 
  
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
    'timeDifference': diff,
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


// Define function that identifies the waterline
/*
declaration:  function identifyWaterFeature(imagem, geometry, scale, bands, selectThreshold(0.0))
    list of parametres:
    • imagem: 
    • geometry:
    • scale:
    • bands: the name of the band in your image
    • selectThreshold: function
    
    global variables: (none)
    libraries needed: (none)
    return value: geometry of the water feature on each image
*/
var histNDWI = [];
var identifyWaterFeature = function (image,geometry,scale,band){

  var internalBandName = 'ndwi';
  var ndwi = image;
  print(ndwi);

  var histogram = ndwi.reduceRegion({
    reducer: ee.Reducer.histogram(),
    geometry: geometry, 
    scale: 10
  });
  histNDWI = histogram.get('ndwi');
  print("histogram",histogram);
  print("histNDWI", histNDWI);
  
  var  threshold = thresholdingAlgorithm(histogram.get('histogram'));
  print ("PRINT 1 - threshold: ", threshold);
  /**
    * Partitions image and reduces to a single vector
  */
  print("PRINT 2 - NDWI", ndwi);
  
  var intertidalFeature = ee.FeatureCollection(ndwi.clip(geometry).lte(threshold)
        .reduceToVectors({ scale: scale, maxPixels: 1e12, eightConnected: false})
        .filter(ee.Filter.eq("label", 1))
        //.limit(1, "count", false) 
        //.first()  
    );
  //print("after reducer", intertidalFeature, typeof(intertidalFeature));
  
  var intertidalFeatureList = intertidalFeature.toList(intertidalFeature.size());
  //print("after to List", intertidalFeatureList, typeof(intertidalFeatureList))
  
  print("PRINT 3 - intertidalFeature", intertidalFeature);
  //
  if(intertidalFeature.getInfo() === null ){
    return null;
  }else{
    return intertidalFeature.geometry();
  }
};

// Define function for gaussianKernel
/*
declaration:  function linearGaussianFilter(size, mean, sigma)
    list of parametres:
    • size: 
    • mean:
    • sigma:

    global variables: (none)
    libraries needed: (none)
    return value: 
*/
var gaussianKernel = function (size, mean, sigma) {
  var gaussianCurve = function (x, mean, sigma)  {
    var divider = ee
      .Number(sigma)
      .multiply(ee.Number(2).multiply(Math.PI).sqrt());
    var exponent = ee.Number(-1).multiply(
      ee
        .Number(x)
        .subtract(mean)
        .pow(2)
        .divide(ee.Number(2).multiply(ee.Number(sigma).pow(2)))
    );

    return ee.Number(1).divide(divider).multiply(exponent.exp());
  };

  var half = ee.Number(size).divide(2).floor();

  var begin = ee.Number(0).subtract(half),
    end = ee.Number(0).add(half);

  // Get the normal distribution Y value for each X
  // in the interval
  var kernel = ee.List.sequence(begin, end).map(function (i) { return gaussianCurve(i, mean, sigma)});

  var sum = kernel.reduce(ee.Reducer.sum());

  // Normalize each value, so that the sum of the list
  // will be equal to one
  var normalizedKernel = kernel.map(function (val) { return ee.Number(val).divide(sum)});

  return normalizedKernel;
};


// Define function that filters with a gaussian curve
/*
declaration:  function linearGaussianFilter(coordinates)
    list of parametres:
    • coordinates: 

    global variables: (none)
    libraries needed: (none)
    return value: 
*/
var linearGaussianFilter = function (coordinates) {
  var samples = 3;
  var mean = 0;
  var sd = 1;
  var coordinateList = ee.List(coordinates);

  // Setup gauss distribution kernel parameters
  var kernelSize = ee.Number(samples);
  var kernelMean = ee.Number(mean);
  var kernelSd = ee.Number(sd);
  var kernel = gaussianKernel(kernelSize, kernelMean, kernelSd);

  var first = coordinateList.reduce(ee.Reducer.first()),
    last = coordinateList.reduce(ee.Reducer.last());

  var sequence = ee.List.sequence(
    ee.Number(0),
    coordinateList.length().subtract(kernelSize)
  );

  var path = sequence.map(function (index) {
    // Take interval of the kernel size to apply the smoothing
    // and zip it to the kernel, so each element in the new list
    // will be a pair of a 2d point and its weight
    var interval = coordinateList
      .slice(ee.Number(index), ee.Number(index).add(kernelSize))
      .zip(kernel);

    // Map the elements, multiplying their axis values by their weight
    var gaussian = interval.map(function (element) {
      // Each element contains a 2d point (0) and a kernel weight (1)
      var asList = ee.List(element);

      var point = ee.List(asList.get(0));
      var weight = ee.Number(asList.get(1));

      // Now we map the two values (each point dimention), multiplying to the weight
      return point.map( function (value) { return ee.Number(value).multiply(weight)});
    });

    // Sum longitude and latitude separately
    var smoothenLong = gaussian
      .map(function (point){ return  ee.List(point).get(0)})
      .reduce(ee.Reducer.sum());
    var smoothenLat = gaussian
      .map(function (point) { return ee.List(point).get(1)})
      .reduce(ee.Reducer.sum());

    // Return final smoothen point
    return ee.List([smoothenLong, smoothenLat]);
  });

  var smoothen = ee.List([]).add(first).cat(path).add(last);

  // return original coordinates if the kernelSize is less than or equal to the length
  // of the given coordinates, otherwise return smoothen coordinates.
  return ee.Algorithms.If(
    coordinateList.size().lte(kernelSize),
    coordinateList,
    smoothen
  );
};


// Calculate time difference and set properties for each feature
function getClosestFeatureImageAndDado(waterlineImage, tideGauge) {
  var closestCollection = tideGauge.map(function(feature) {
    return calculateTimeDifferenceAndSetProperties(waterlineImage, feature);
  });

  var sortedCollection = closestCollection.sort('timeDifference');
  var closestFeatureImage = sortedCollection.first();
  var tideValue = closestFeatureImage.get('Dado');
  var dateImage = closestFeatureImage.get('dateImage');
  var dateGauge = closestFeatureImage.get('dateGauge');
  
  return {
    closestFeatureImage: closestFeatureImage,
    tideValue: tideValue,
    dateImage: dateImage,
    dateGauge: dateGauge
  };
}


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

var evaluated = ee.FeatureCollection([]);
//var SDB_intertidal = function(pathImage, tg, savePath) {
var assetList = ee.data.listImages(pathImage).images;
var idsList = assetList.map(function(img) { return img.id}); 
var imageCollection = ee.ImageCollection(idsList);
print("Before datetime:", ee.Date(imageCollection.first().get("system:time_start")).format(null, "UTC"));
var imageCollectionWithDatetime = imageCollection.map(addDatetimeProperty);
print('Image Collection with Datetime:', imageCollectionWithDatetime);
print("After datetime:", (imageCollectionWithDatetime.first().get("datetime"))); 

var numImages = imageCollectionWithDatetime.size();
//var mergedPointsWithTide = ee.FeatureCollection([]);
var closestFeatureImage; 

print('-- begin waterline print --');


//////////////////// FOR LEGEND
// Create a legend panel
var legend = ui.Panel({
style: {
  position: 'bottom-right',
  padding: '8px 15px'
}
});

// Create a title for the legend
var legendTitle = ui.Label({
value: 'Waterline Legend',
style: {fontWeight: 'bold', fontSize: '18px', margin: '0 0 10px 0', padding: '0'}
});
legend.add(legendTitle);

// Define a function to add a color and label to the legend
var addColorAndLabel = function(color, label) {
// Create a label
var colorBox = ui.Label({
  style: {
    backgroundColor: color,
    padding: '8px',
    margin: '0 0 4px 0'
  }
});

// Create a panel with the color box and the label
var legendEntry = ui.Panel({
  widgets: [colorBox, ui.Label(label)],
  layout: ui.Panel.Layout.Flow('horizontal'),
  style: {margin: '0 0 4px 0'}
});

legend.add(legendEntry);
};

///////////////////// END LEGEND


for (var i = 0; i < numImages.getInfo(); i++) {
  
  var ndwi = ee.Image(imageCollectionWithDatetime.toList(numImages).get(i));
  Map.addLayer(ndwi, {min:-1, max:1}, "NDWI" +i, false);
  var datetime = ee.Date(ndwi.get("datetime")).format().getInfo();  // Get datetime as string
  
  var result = getClosestFeatureImageAndDado(ndwi, tideGauge);
  print(result)
  var closestFeatureImage = result.closestFeatureImage;
  var tideValue = result.tideValue;
  var dateImage = result.dateImage;
  var dateGauge = result.dateGauge;
  
  print("PRINT 3.1 - tide value closest to the image date", tideValue)
      
  var hist_NDWI1 = ui.Chart.image.histogram({image:ndwi.select('ndwi'), scale: 10, maxPixels: 10e9})
    .setOptions({title: 'Histograma  NWI sem corte '+i});
  print(hist_NDWI1);
  
  var band = "ndwi";
  var waterSegment = identifyWaterFeature(ndwi, geometry, scale, band); //, selectThreshold(-1.0)
  print("PRINT 4 - waterSegment: "+i, waterSegment);
  
  if(waterSegment !== null){
    Map.addLayer(waterSegment,
                 {'color': 'black'},
                 'watersegment: '+i, false);
                 
    var polygons = waterSegment.coordinates();
    print("polygons", polygons);
    
    // Function to filter polygons with more than 25 vertices
    var filteredPolygons = polygons.map(function(polygon) {
      var numVertices = ee.List(ee.List(polygon).get(0)).length();
      return ee.Algorithms.If(numVertices.gt(20), polygon, null);
    }).removeAll([null]);
    
    waterSegment = ee.Geometry.MultiPolygon(filteredPolygons);
    print("PRINT 4.1 - waterSegment: "+i, waterSegment);
    
    polygons = waterSegment.coordinates();
    var processedPolygons = polygons.map(function(polygon) {
      var size = ee.List(polygon).size();
      return ee.Algorithms.If(size.gt(1), [ee.List(polygon).get(0)], polygon);
    });

    // Create a new MultiPolygon with the processed polygons
    var waterSegment = ee.Geometry.MultiPolygon(processedPolygons); 
    print("PRINT 4.2 - waterSegment: "+i, waterSegment);
    Map.addLayer(ee.FeatureCollection(waterSegment), {color: 'red'}, 'Filtered waterSegment: '+i, false);

    var waterline = ee.FeatureCollection([]);

    for (var ii=0; ii < waterSegment.coordinates().size().getInfo(); ii++) {
      
      var coordinates = waterSegment.coordinates().get(ii);
      var wl = ee.FeatureCollection(ee.Geometry.MultiLineString(coordinates)
        .intersection(geometry.buffer(-10)));
        
      wl = ee.FeatureCollection(ee.Geometry.MultiLineString(
        linearGaussianFilter(wl.first().geometry().coordinates())));
      
      var updatedWl = wl.map(function(feature) {
        return ee.Feature(feature.geometry(), feature.toDictionary())
          .set('id', 'WL_' + datetime + '_' + ii)
          .set('dateTime', datetime)
          .set('tideValue', tideValue)
          .set('dateGauge', ee.Date(dateGauge).format('yyyy-MM-dd'))
          .set('dateImage', ee.Date(dateImage).format('yyyy-MM-dd')); 
      });
      
      waterline = waterline.merge(updatedWl);
       
      }
      
    print("PRINT 5 - Waterline "+i, waterline);
    
    
    var convertToLineString = function(feature) {
      var coords = ee.List(ee.List(feature.geometry().coordinates()).get(0));
      var proj = 'EPSG:4326';
      return ee.Feature(ee.Geometry.LineString(coords, proj), feature.toDictionary()); //{"system:time_start": ee.Date(ndwi.get("datetime")).format(null)}
    };
    
    var waterlines = waterline.map(convertToLineString);
    
    print ("PRINT 6 - waterlines: " + i, waterlines);
    
    //evaluated.push(waterlines);
    
    var colors = ['#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#00FFFF', '#FF00FF', '#FFA500', '#800080', '#008080', '#000000'];
    var color = colors[i % colors.length];
    
    Map.addLayer(waterline,
                  {'color': color }, //{'color': 'blue'},
                 'waterline: ' + i);
                 
    print("-- end waterline print --");
    
    // Add color and datetime to the legend
    addColorAndLabel(color, datetime);
  
  print(noncense)  
  }else{
    print("Error image: "+ i);
  }
  evaluated = evaluated.merge(waterlines);
  
  /*
  var id = ee.Date(dateImage).format('yyyyMMddHHmmss').getInfo();
  //to export each
  Export.table.toDrive({
    collection: waterlines,
    description: 'waterlines_Delta' + id,
    fileFormat: 'SHP',
    folder:'Intertidal_Zones'
  });
  */
}

// Add the legend to the map
Map.add(legend);

//to export the whole collection
Export.table.toDrive({
  collection: evaluated,
  description: 'waterlines_Delta',
  fileFormat: 'SHP',
  folder:'Intertidal_Zones'
});
