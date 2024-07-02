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

var pathImage = 'projects/ee-index-images/assets/amor'; // change this to the path of the asset containing the NWI data
var pathTideg = 'projects/ee-tide-gauge/assets/'; // change this to the path of the asset containing the tide gauge data
var nameTideg = 'tideVale'; // name of the tide gauge observations file
var savePath = 'projects/ee-tide-gauge/assets/'; // change this to the path of the asset where you want to save the results
var tideGauge = ee.FeatureCollection(pathTideg + nameTideg); // loading tide gauge observations as a ee.FeatureCollection
print('First 10 elements of tide gauge:',tideGauge.limit(10));

// algumas definições
var scale = 10;
var geometry = ee.Geometry.Polygon(
[[[-44.35859636626557,-2.5275481643998887],
              [-44.35859636626557,-2.533703846828047],
              [-44.34756866598747,-2.533703846828047],
              [-44.34756866598747,-2.5275481643998887],
              [-44.35859636626557,-2.5275481643998887]]]);
Map.centerObject(geometry, 15);
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
  var dateTime = ee.Date.parse('dd/MM/yyyy HH:mm:ss', dateTimeStr); //yyyy-MM-dd HH:mm:ss 
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
  var datetimeStr = id.slice(37, 52);  // Extract the datetime string from the ID 
  
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
  //print(ndwi);

  var histogram = ndwi.reduceRegion({
    reducer: ee.Reducer.histogram(),
    geometry: geometry, 
    scale: 10
  });
  histNDWI = histogram.get('ndwi');
  //print("histogram",histogram);
  //print("histNDWI", histNDWI);
  
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


// Masking filter - ndwi
// Input: image collection with NDWI band added 
// Output: image collection with land areas masked
// Description: Mask all land pixels in image NDWI 
//  Values correspond to the following ranges:
//    0.2 – 1 – Water surface;
//    0.0 – 0.2 – Flooding, humidity;
//    -0.3 – 0.0 – Moderate drought, surfaces without water;
//    -1 – -0.3 – Dry, water-free surfaces.

var waterLandFilter = function (image){
  return image.updateMask(image.lt(0.2))}; 


// Function to filter polygons with more than 20 vertices
// Input: multiPolygons 
// Output: multiPolygons
// Description: Mask all polygons from a feature collection of polygons that have less then 20 vertices
var filterPolygons = function(polygon) {
  var numVertices = ee.List(ee.List(polygon).get(0)).length();
  return ee.Algorithms.If(numVertices.gt(20), polygon, null);
};

// Function to get the polygon with most vertices
// Input: multiPolygons 
// Output: Polygon
// Description: Get the biggest polygon from a freature of multipolygons
var getPolygonWithMostVertices = function(feature) {
  var polygons = feature.geometry().geometries();
  
  var polygonsWithVertexCount = polygons.map(function(polygon) {
    var numVertices = ee.Geometry(polygon).coordinates().flatten().length();
    return ee.Feature(ee.Geometry(polygon)).set('numVertices', numVertices);
  });
  var polygonsCollection = ee.FeatureCollection(polygonsWithVertexCount);
  var sortedPolygons = polygonsCollection.sort('numVertices', false);
  var largestPolygon = sortedPolygons.first().geometry();
  
  return ee.Feature(largestPolygon);
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

var evaluated = ee.FeatureCollection([]);
//var SDB_intertidal = function(pathImage, tg, savePath) {
var assetList = ee.data.listImages(pathImage).images;
var idsList = assetList.map(function(img) { return img.id}); 
var imageCollection = ee.ImageCollection(idsList);
//print("Before datetime:", ee.Date(imageCollection.first().get("system:time_start")).format(null, "UTC"));
var imageCollectionWithDatetime = imageCollection.map(addDatetimeProperty);
print('Image Collection with Datetime:', imageCollectionWithDatetime);
//print("After datetime:", (imageCollectionWithDatetime.first().get("datetime"))); 

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
  var useCase = i;
  switch(useCase){
  case 0: 
      var toIntercept = ee.Geometry.Polygon(
        [[[-44.354473095265185, -2.533902608682677],
          [-44.35028884920195, -2.5329593939091053],
          [-44.347231130924975, -2.5306442274622887],
          [-44.34830401453093, -2.528992351535467],
          [-44.34862587961272, -2.5293353398551295],
          [-44.34840849151947, -2.5299765222455926],
          [-44.34939554443695, -2.530973331041826],
          [-44.35173443069794, -2.5318308003607646],
          [-44.35312917938568, -2.531948702347744],
          [-44.35418060531952, -2.5320558859629627],
          [-44.3547170471225, -2.53215235120908],
          [-44.355532438663026, -2.5322916898852315],
          [-44.356562406924745, -2.5324095918302785],
          [-44.3566804241214, -2.5339423161391683]]]);
  break;
  case 1:
    var toIntercept = ee.Geometry.Polygon(
      [[[-44.35670958562146, -2.532389655284352],
            [-44.35639844937573, -2.5327647977711476],
            [-44.354853496983154, -2.5327647977711476],
            [-44.34942470593701, -2.531993075965908],
            [-44.34724675221692, -2.5299780224189514],
            [-44.348308906986816, -2.5292920460288295],
            [-44.348555670216186, -2.5296277535035903],
            [-44.34827672047864, -2.530217264333165],
            [-44.34874878926526, -2.530742464664666],
            [-44.34945689244519, -2.5312033545762067],
            [-44.35016499562512, -2.5315356239456417],
            [-44.35122715039502, -2.531685681052305],
            [-44.35213910146008, -2.531749991235556],
            [-44.353083239033325, -2.5318786115924765],
            [-44.354563818409545, -2.5321572889886843],
            [-44.35551868481885, -2.532264472586652],
            [-44.35594783826123, -2.5322430358677606],
            [-44.35616241498242, -2.532489558113509],
            [-44.35628043217908, -2.5324145296089036],
            [-44.35637699170361, -2.532200162428935]]]);
  break;
  case 2:
    var toIntercept = ee.Geometry.Polygon(
      [[[-44.357508644626684, -2.533890312691725],
            [-44.354729876087255, -2.5338795943457777],
            [-44.34941397710095, -2.5321290893497235],
            [-44.347053633167846, -2.530328403539231],
            [-44.348609314396484, -2.528206163480651],
            [-44.34940324826489, -2.5282383186589596],
            [-44.34942470593701, -2.5296671894121343],
            [-44.35463892026196, -2.5318858923454983],
            [-44.35617314381848, -2.532282471696753],
            [-44.35640917821179, -2.531993075965908]]]);
  break;
  case 3:
    var toIntercept = ee.Geometry.Polygon(
      [[[-44.357508644626684, -2.533890312691725],
            [-44.354729876087255, -2.5338795943457777],
            [-44.34941397710095, -2.5321290893497235],
            [-44.347053633167846, -2.530328403539231],
            [-44.348609314396484, -2.528206163480651],
            [-44.34940324826489, -2.5282383186589596],
            [-44.34942470593701, -2.5296671894121343],
            [-44.35463892026196, -2.5318858923454983],
            [-44.35617314381848, -2.532282471696753],
            [-44.35640917821179, -2.531993075965908]]]);
  break;    
  case 4:
    var toIntercept = ee.Geometry.Polygon(
        [[[-44.35115717367942, -2.532539800404368],
          [-44.34694074110801, -2.5307498335364564],
          [-44.34882901625449, -2.528413226361354],
          [-44.34931181387717, -2.5286811861134613],
          [-44.35035251097495, -2.5299673921530843],
          [-44.35162924246604, -2.53051402933365],
          [-44.35362480597312, -2.5313822173231304],
          [-44.35518048720176, -2.5322396863715055],
          [-44.355706200168676, -2.5323575883212723],
          [-44.35560964064414, -2.53248620861789]]]);
  break;
  case 5:
    var toIntercept = ee.Geometry.Polygon(
        [[[-44.35327947586651, -2.5328621052516422],
          [-44.35170770138378, -2.53299608467481],
          [-44.34876250164798, -2.5318459512126874],
          [-44.34733020203403, -2.5301685263465825],
          [-44.34808658497623, -2.5295468603497806],
          [-44.34855865376285, -2.530490077607302],
          [-44.3492882146149, -2.531052792121629],
          [-44.3504040135651, -2.531545837019105],
          [-44.35127841370395, -2.5317494859437595],
          [-44.35235666172794, -2.532044240909599],
          [-44.35275899308017, -2.5320335225483723],
          [-44.35289310353092, -2.5320335225483723],
          [-44.35306476490787, -2.532081755173174],
          [-44.35326324837497, -2.5323121999112606]]]);
//Map.addLayer(toIntercept, {}, "toIntercept");
}  

  var ndwi = ee.Image(imageCollectionWithDatetime.toList(numImages).get(i));
  Map.addLayer(ndwi, {min:-1, max:1}, "NDWI" +i, false);
  
  ///// CREATE POLYGON OF "NON WATER FEATURE" TO INTERSEPT:
  //var waterLandImage = ndwi.map(waterLandFilter); 
  
  
  var histogram = ndwi.reduceRegion({
    reducer: ee.Reducer.histogram(),
    geometry: geometry, 
    scale: 10
  });
  histNDWI = histogram.get('ndwi');
  var threshold = thresholdingAlgorithm(histogram.get('histogram'));
  var minPixelValue = ndwi.reduceRegion({reducer: ee.Reducer.min(),geometry: geometry,scale: 10, maxPixels: 1e9});
  var maxPixelValue = ndwi.reduceRegion({reducer: ee.Reducer.max(),geometry: geometry,scale: 10, maxPixels: 1e9});
  var stdDevPixelValue = ndwi.reduceRegion({reducer: ee.Reducer.stdDev(),geometry: geometry,scale: 10, maxPixels: 1e9});
  stdDevPixelValue = stdDevPixelValue.get('ndwi').getInfo();
  print("stdDev value", stdDevPixelValue);  
  var THajustment = ee.Number(stdDevPixelValue).add(1);
  var thresholdValue = ee.Number(threshold).subtract(ee.Number(THajustment).multiply(ee.Number(threshold))); 
  //print(thresholdValue);

  var waterLandImageFeature = ee.FeatureCollection(ndwi.clip(geometry).lt(thresholdValue)
      .reduceToVectors({ scale: scale, maxPixels: 1e12, eightConnected: false})
      .filter(ee.Filter.eq("label", 1))
  );
  print("waterLandImageFeature", waterLandImageFeature);

  var filteredWaterLandImageFeature = ee.Geometry.Polygon(ee.List(getPolygonWithMostVertices(waterLandImageFeature)
    .geometry().coordinates()).get(0));
  Map.addLayer(filteredWaterLandImageFeature, {'color': 'grey'}, 'waterLandImageFeature: '+i, false);

  var datetime = ee.Date(ndwi.get("datetime")).format().getInfo();  // Get datetime as string
  
  var result = getClosestFeatureImageAndDado(ndwi, tideGauge);
  var closestFeatureImage = result.closestFeatureImage;
  var tideValue = result.tideValue;
  var dateImage = result.dateImage;
  var dateGauge = result.dateGauge;
  
  var hist_NDWI1 = ui.Chart.image.histogram({image:ndwi.select('ndwi'), scale: 10, maxPixels: 10e9})
    .setOptions({title: 'Histograma  NWI sem corte '+i});
  print(hist_NDWI1);
  
  var band = "ndwi";
  var waterSegment = identifyWaterFeature(ndwi, geometry, scale, band); //, selectThreshold(-1.0)
  //print("PRINT 4 - waterSegment: "+i, waterSegment);
  
  if(waterSegment !== null){
    Map.addLayer(waterSegment,
                 {'color': 'black'},
                 'watersegment: '+i, false);
                 
    var polygons = waterSegment.coordinates();
    print("polygons", polygons);
    
    var filteredPolygons = polygons.map(filterPolygons).removeAll([null]);
    
    waterSegment = ee.Geometry.MultiPolygon(filteredPolygons);
    //print("PRINT 4 - waterSegment: "+i, waterSegment);
    
    polygons = waterSegment.coordinates();
    var processedPolygons = polygons.map(function(polygon) {
      var size = ee.List(polygon).size();
      return ee.Algorithms.If(size.gt(1), [ee.List(polygon).get(0)], polygon);
    });

    // Create a new MultiPolygon with the processed polygons
    var waterSegment = ee.Geometry.MultiPolygon(processedPolygons); 
    print("PRINT 4 - waterSegment: "+i, waterSegment);  
    Map.addLayer(ee.FeatureCollection(waterSegment), {color: 'red'}, 'Filtered waterSegment: '+i, false);
    //print(nonsense)

    var waterline = ee.FeatureCollection([]);

    for (var ii=0; ii < waterSegment.coordinates().size().getInfo(); ii++) {
      
      var coordinates = waterSegment.coordinates().get(ii);
      //print("PRINT 4.1 - waterline before difference", coordinates)
      var wl = ee.FeatureCollection(ee.Geometry.MultiLineString(coordinates)
        .difference(toIntercept)); //intersection //.buffer(-5)
      //print("PRINT 4.2 - waterline after difference: "+ii, wl);
      //print(ee.List(wl.first()))
      
      var checkCoordinates = wl.first().geometry().coordinates();
      
      if (ee.List(checkCoordinates).size().gt(0).getInfo()) {
        var type2test = ee.Algorithms.ObjectType(ee.List(ee.List(checkCoordinates).get(0)).get(0));
          if (type2test.equals('List').getInfo()) {
            var coordsList = ee.List(ee.List(checkCoordinates)).get(-1);
            wl = ee.FeatureCollection(ee.Geometry.MultiLineString(coordsList));
            print("PRINT 4.4 - waterline after correction: ", wl);  // Print the corrected wl for debugging
          } else {
            print("Error: type2test is not 'List', it is", type2test);
        } 
        
      wl = ee.FeatureCollection(ee.Geometry.MultiLineString(
        linearGaussianFilter((wl).geometry().coordinates()))); 
      print("PRINT 4.4 - waterline after smooth: "+ii, wl);
      Map.addLayer(wl, {}, "tst fdp", false);
      
      var updatedWl = wl.map(function(feature) {
        return ee.Feature(feature.geometry(), feature.toDictionary())
          .set('id', 'WL_' + datetime + '_' + ii)
          .set('dateTime', datetime)
          .set('tideValue', tideValue)
          .set('dateGauge', ee.Date(dateGauge).format('yyyy-MM-dd'))
          .set('dateImage', ee.Date(dateImage).format('yyyy-MM-dd')); 
      });
      
      waterline = waterline.merge(updatedWl);
      
      } else {
        print("Error: coordinates list is empty.");
      } 
    }
      //print("PRINT 5 - Waterline "+i, waterline);
    
    
    var convertToLineString = function(feature) {
      var coords = ee.List(ee.List(feature.geometry().coordinates()).get(0));
      var proj = 'EPSG:4326';
      return ee.Feature(ee.Geometry.LineString(coords, proj), feature.toDictionary()); //{"system:time_start": ee.Date(ndwi.get("datetime")).format(null)}
    };
    
    var waterlines = waterline.map(convertToLineString);
    
    //print ("PRINT 6 - waterlines: " + i, waterlines);
    
    var colors = ['#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#00FFFF', '#FF00FF'];
    //'#FFA500', '#800080', '#008080', '#000000'];
    var color = colors[i % colors.length];
    
    Map.addLayer(waterline,{'color': color },'waterline: ' + i, false);
    print("-- end waterline print --");
    
    addColorAndLabel(color, datetime); // Add color and datetime to the legend
  }else{
    print("Error image: "+ i);
  }
  evaluated = evaluated.merge(waterlines);
}

// Add the legend to the map
Map.add(legend);

//to export the whole collection
Export.table.toDrive({
  collection: evaluated,
  description: 'waterlines_Amor',
  fileFormat: 'SHP',
  folder:'Intertidal_Zones'
});
