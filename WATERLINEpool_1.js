var index = 'ndwi'; // ndwi, mndwi_sharp, awei
var thLand = 0; // -0.2 para macromaré, 0 ou 0.2 para os demais
var useCase = 1; //1= tauranga, 2 = amor, 3 = vale, 4 = bsul, 5 = delta, 6 = aterro

//use case selection
var geometry = null;
switch(useCase){
  case 1:
    //Entradas para TAURANGA
    var img0 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20190116T221559_20190116T221559_T60HVD');
    var img1 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20190205T221559_20190205T221600_T60HVD');
    var img2 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20190215T221309_20190215T221543_T60HVD');
    var img3 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20190225T221559_20190225T221559_T60HVD');
    var img4 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20190302T221601_20190302T221555_T60HVD');
    var img5 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20190406T221609_20190406T221604_T60HVD');
    var img6 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20190501T221601_20190501T221604_T60HVD');
    
    var imageCollection = ee.ImageCollection.fromImages([img0, img1, img2, img3, img4, img5, img6]);  
    
    var geometry = ee.FeatureCollection('projects/ee-tide-gauge/assets/tauranga');
    Map.centerObject(geometry, 11);
    
  break;
  
  case 2:
    //Entradas para PRAIA DO AMOR
    var img0 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20230718T132239_20230718T132238_T23MNT');
    var img1 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20230723T132241_20230723T132237_T23MNT'); //boa só Praia do Amor
    var img2 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20230812T132241_20230812T132237_T23MNT'); //boa só Praia do Amor
    var img3 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20230827T132239_20230827T132236_T23MNT');
    var img4 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20231026T132239_20231026T132233_T23MNT');
    var img5 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20231031T132231_20231031T132235_T23MNT');
    
    var imageCollection = ee.ImageCollection.fromImages([img0,img1, img2, img3, img4, img5]);  
    
    var geometry = ee.Geometry.Polygon(
    [[[-44.35859636626557,-2.5275481643998887],
              [-44.35859636626557,-2.533703846828047],
              [-44.34756866598747,-2.533703846828047],
              [-44.34756866598747,-2.5275481643998887],
              [-44.35859636626557,-2.5275481643998887]]]);
    Map.centerObject(geometry, 15);

  break;
  
  case 3:
    //Entradas para PRAIA DA VALE
    var img0 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20230718T132239_20230718T132238_T23MNT');
    var img1 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20230827T132239_20230827T132236_T23MNT');
    var img2 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20231026T132239_20231026T132233_T23MNT');
    var img3 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20231031T132231_20231031T132235_T23MNT');
    
    var imageCollection = ee.ImageCollection.fromImages([img0,img1, img2, img3]);  
    
    var geometry = ee.Geometry.Polygon(
    [[[-44.37871414728457,-2.5585251544134593],
              [-44.368659677502194,-2.564403867621195],
              [-44.362658172283204,-2.560510683083635],
              [-44.37520677410478,-2.55373651379837],
              [-44.37871414728457,-2.5585251544134593]]]);
    Map.centerObject(geometry, 15);   
  break;
  case 4:
    //Entradas para BABITONGA / BARRA DO SUL
    var img0 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20230320T132239_20230320T132815_T22JGR');
    var img1 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20230330T132239_20230330T132858_T22JGR');
    var img2 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20230404T132231_20230404T132232_T22JGR');
    var img3 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20230419T132239_20230419T132235_T22JGR');
    var img4 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20230514T132231_20230514T132233_T22JGR');
    var img5 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20230519T132239_20230519T132237_T22JGR');
    var img6 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20230603T132241_20230603T132236_T22JGR');
    var img7 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20230608T132239_20230608T132237_T22JGR');
    var img8 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20230628T132239_20230628T132237_T22JGR');
    var img9 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20230713T132241_20230713T132606_T22JGR');
    var img10 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20230807T132239_20230807T132237_T22JGR');
    
    var imageCollection = ee.ImageCollection.fromImages([img0,img1, img2, img3, img4, img5, img6, img7, img8, img9, img10]);  
    
    var geometry = ee.Geometry.Polygon(
        [[[-48.62734331348832, -26.44521206545918],
          [-48.61755861500199, -26.45381887874305],
          [-48.6127520964473, -26.455739954609705],
          [-48.60974802235062, -26.4591978103976],
          [-48.60485567310746, -26.458429406975334],
          [-48.59558595875199, -26.455125213819553],
          [-48.59790338734086, -26.449054472293565],
          [-48.60202326038773, -26.45266621783879],
          [-48.61052049854691, -26.448209153787833],
          [-48.61300253233181, -26.443604534281864],
          [-48.61746572813259, -26.440761031028728]]]);
    Map.centerObject(geometry, 14); 
  
  break;
  case 5:
    //Entradas para BABITONGA / DELTA DE ENCHENTE
    var img0 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20230320T132239_20230320T132815_T22JGS');
    var img1 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20230330T132239_20230330T132858_T22JGS');
    var img2 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20230404T132231_20230404T132232_T22JGS');
    var img3 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20230419T132239_20230419T132235_T22JGS');
    var img4 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20230514T132231_20230514T132233_T22JGS');
    var img5 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20230519T132239_20230519T132237_T22JGS');
    var img6 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20230603T132241_20230603T132236_T22JGS');
    var img7 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20230608T132239_20230608T132237_T22JGS');
    var img8 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20230628T132239_20230628T132237_T22JGS');
    var img9 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20230713T132241_20230713T132606_T22JGS');
    var img10 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20230807T132239_20230807T132237_T22JGS');
    
    var imageCollection = ee.ImageCollection.fromImages([img0,img1, img2, img3, img4, img5, img6, img7, img8, img9, img10]);  
    
    var geometry = ee.Geometry.Polygon(
    [[[-48.6064338684082,-26.22352298843246],
                [-48.55304718017578,-26.22352298843246],
                [-48.55304718017578,-26.185018250078308],
                [-48.6064338684082,-26.185018250078308],
                [-48.6064338684082,-26.22352298843246]]]);
    Map.centerObject(geometry, 14); 
    
  break;
  case 6:
    //Entradas para BABITONGA / ATERRO
    var img0 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20230320T132239_20230320T132815_T22JGR');
    var img1 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20230330T132239_20230330T132858_T22JGR');
    var img2 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20230404T132231_20230404T132232_T22JGR');
    var img3 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20230419T132239_20230419T132235_T22JGR');
    var img4 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20230514T132231_20230514T132233_T22JGR');
    var img5 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20230519T132239_20230519T132237_T22JGR');
    var img6 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20230603T132241_20230603T132236_T22JGR');
    var img7 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20230608T132239_20230608T132237_T22JGR');
    var img8 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20230628T132239_20230628T132237_T22JGR');
    var img9 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20230713T132241_20230713T132606_T22JGR');
    var img10 = ee.Image('COPERNICUS/S2_SR_HARMONIZED/20230807T132239_20230807T132237_T22JGR');
    
    var imageCollection = ee.ImageCollection.fromImages([img0,img1, img2, img3, img4, img5, img6, img7, img8, img9, img10]);  
    
    var geometry = ee.Geometry.Polygon(
    [[[-48.74015808105469,-26.37372339295363],
                [-48.65947723388672,-26.37372339295363],
                [-48.65947723388672,-26.31895972813747],
                [-48.74015808105469,-26.31895972813747],
                [-48.74015808105469,-26.37372339295363]]]);
    Map.centerObject(geometry, 13);
}

print('Image Collection:', imageCollection);

// Mask out land and clouds
function mask(image) {
  var qa = image.select('QA60');

  // Bits 10 and 11 are clouds and cirrus, respectively.
  var cloudBitMask = 1 << 10; //var cloudBitMask = ee.Number(2).pow(10).int();
  var cirrusBitMask = 1 << 11; //var cirrusBitMask = ee.Number(2).pow(11).int();

  // Both flags should be set to zero, indicating clear conditions.
  var mask_ = qa.bitwiseAnd(cloudBitMask).eq(0)
    .and(qa.bitwiseAnd(cirrusBitMask).eq(0));

  return image.updateMask(mask_).divide(10000);
}

var Bay_Plenty_collection = imageCollection;
//Bay_Plenty_collection = Bay_Plenty_collection.map(mask);
print ('initial collection', Bay_Plenty_collection.size().getInfo());

function reflec_corr (image){
  var opticalBands = image.select('B.*').multiply(0.0001); //applying the scale factor of Sentinel-2 collection
    return image
    .addBands(opticalBands, null, true);
  }
var Bay_Plenty_collection = Bay_Plenty_collection.map(reflec_corr);
print ('scaled collection',Bay_Plenty_collection);

var vis_rgb = {
  max: 0.25,
  bands: ['B4', 'B3', 'B2'], // to visualize in true color
  };
var singleBandVis = {
              'min': -0.5,
              'max': 1,
              };
Map.addLayer(Bay_Plenty_collection.first(), vis_rgb, "first image masked");

///// PANSHARPENING
var geeSharp = require('users/aazuspan/geeSharp:geeSharp'); // Import the geeSharp module
// After analysing the charts, choose the band that showed the bigest r^2
function sharpened (image) {
  var sharp1 = geeSharp.sharpen(image.select(['B11']), image.select(['B8'])).rename ('B11_sharp');
  var sharp2 = geeSharp.sharpen(image.select(['B12']), image.select(['B8'])).rename ('B12_sharp');
  
  return image.addBands([sharp1,sharp2])}

var Bay_Plenty_collection = Bay_Plenty_collection.map(sharpened);
print('sharpened collection',Bay_Plenty_collection);

//// indices
var indexFunction = function(image) { 
 var ndwi = image.normalizedDifference(['B3', 'B8']).rename ('ndwi'); //Mc Feeters, 1996
 var mndwi_sharp = image.normalizedDifference(['B3', 'B11_sharp']).rename ('mndwi_sharp');
 var mndwi = image.normalizedDifference(['B3', 'B11']).rename ('mndwi');
 //var ndvi = image.normalizedDifference(['B3', 'B8']).rename ('ndwi');
 var awei = image.expression('(B+ (2.5*G) -1.5*(N+S1) -(0.25*S2)) ',{ //Feyisa etal, 2014 (4 * (G - S)) - (0.25 * N + 2.75 * S)
   B: image.select('B2'),
   G: image.select('B3'), 
   S1: image.select('B11_sharp'), 
   S2: image.select('B12_sharp'),
   N: image.select('B8'),
   }).rename('awei');
   
  return image.addBands([ndwi, mndwi_sharp,mndwi,awei]);
};
var clip_image = function (image){
  return image.clip(geometry)}; // corta para a área de interesse
  
var mask_land = function (image){
  var ndwi = image.select('mndwi_sharp');
  return image.updateMask(ndwi.gte(thLand))}; // deixa apenas valores maiores ou iguais a 0.2 //-0,2 para macromaré

///// APLICANDO AS FUNÇÕES E CORTANDO NA GEOMETRIA
var NWI = Bay_Plenty_collection.map(indexFunction);// calcula o NDWI;  
print('pós NDWI', NWI);
var NWINoMask=NWI.map(clip_image); // Faz o corte

//// MÁSCARA DE ÁGUA/TERRA
var NWI=NWINoMask.map(mask_land);  //// aplica máscara a aprtir do valor do histograma
Map.addLayer(NWI.first().select('mndwi_sharp'),singleBandVis, "NDWI-Masked");
Map.addLayer(NWINoMask.first().select('mndwi_sharp'),singleBandVis, "NDWI-No Masked",false);

var aweiVis = {
   'max': 0.2,
   'min': -0.2,
};
Map.addLayer(NWI.first().select(index),aweiVis, "AWEI index");

// Visualizar AWEI histograma
var hist_NDWI_Mask = ui.Chart.image.histogram({image:NWI.first().select('mndwi_sharp'), region: geometry, scale: 11})
  .setOptions({title: 'Histograma mNDWI com máscara'});
print (hist_NDWI_Mask);

var palette = ['blue','yellow'];
var vis_params = {
              'min': -1,
              'max': 1,
              'dimensions': 500,
              'palette':palette,             
              };

var hist_AWEI_Mask = ui.Chart.image.histogram({image:NWI.select(index).first(), region: geometry, scale: 11})
  .setOptions({title: 'Histograma AWEI com máscara - ' +index});
print (hist_AWEI_Mask);

var NWI_STD = NWI.select(index).reduce(ee.Reducer.stdDev()); // Agora só tem uma imagem que mostra o STD dos NDWI
Map.addLayer(NWI_STD, singleBandVis,'STD image - '+index);

var hist_NWI_STD = ui.Chart.image.histogram({image:NWI_STD, region: geometry, scale: 11})
  .setOptions({title: 'Histograma AWEI STD'});
print (hist_NWI_STD);

// UTILIZANDO A METODOLOGIA DE THRESHOLD OTSU:
var histogram = NWI_STD.reduceRegion({
  reducer: ee.Reducer.histogram(),
  geometry: geometry, 
  scale: 10
});

if (index === 'ndwi') {
    var stdDev = 'ndwi_stdDev';
} else { if (index === 'mndwi') {
    var stdDev = 'mndwi_stdDev';
    } else {
      var stdDev = 'awei_stdDev';
    }
}
var AWEI_stdDev = histogram.get(stdDev);
print(AWEI_stdDev);

var otsu = function(histogram) {
  var counts = ee.Array(ee.Dictionary(AWEI_stdDev).get('histogram'));
  var means = ee.Array(ee.Dictionary(AWEI_stdDev).get('bucketMeans'));
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

  print(ui.Chart.array.values(ee.Array(bss), 0, means));

  // Return the mean value corresponding to the maximum BSS.
  return means.sort(bss).get([-1]);
};

var threshold = otsu(histogram.get('histogram'));
print('the threshold is:    ', threshold);

///// Máscara a partir do histograma do desvio pasdrão
// valores maiores que X*STD (retirar X do histograma do NWI_STD)
print("NWI_STD", NWI_STD);
var stdMasked = NWI_STD.updateMask(NWI_STD.gte(threshold)); //EQ("="), GTE(">="), GT(">"), LT("<"), LTE("<="); // Threshold para Barra do Sul: 0.05
var zones=NWI_STD.gte(threshold); //threshold is mean value from the histogram
var zones = zones.updateMask(zones.neq(0));
//print("zones...", zones);
Map.addLayer(zones, vis_params, 'stdMasked', false); 
// print("zones", zones)
// Visualizar histograma pós corte
var hist_intermare = stdMasked.reduceRegion({
  reducer: ee.Reducer.histogram(),
  geometry: geometry, 
  scale: 10
});

var hist_intermare = ui.Chart.image.histogram({image:stdMasked, region:geometry, scale:11})
  .setOptions({title: 'Histograma intermaré'});
print (hist_intermare);

///// Transformar para vetores
var vectors = zones.addBands(NWI_STD).reduceToVectors({ 
  crs: 'EPSG:4326',
  //crsTransform: [60, 0, 399960, 0, -60, 5900020],
  scale: 10,
  geometryType:'polygon',
  labelProperty: 'stdMasked',
  eightConnected: false,
  geometry: geometry,
  maxPixels: 100e9,
  geometryInNativeProjection: true,
  reducer: ee.Reducer.mean()
});
//print("número de vetores identificados: ", (vectors.getInfo())); 

///// Função para cortar de acordo com os vetores
var clip_image2 = function(image){
  return image.clip(vectors);
};

var NWI2 = Bay_Plenty_collection.map(indexFunction).select(index);

print("Clip das imagens", NWI2);
var intertidal_zones = NWI2.map(clip_image2);

// Displaying
var palette = ['white','blue'];
var vis_params = {
              'min': -0.5,
              'max': 1,
              'dimensions': 500,
              'palette':palette,             
              };
Map.addLayer(intertidal_zones.select(index), vis_params, 'intertidal zone '+index);

/* Exportar vetores
var task = Export.table.toDrive({
  collection: vectors,
  description:'Intertidal_Zones_Linguado_AWEI',
  folder: 'Intertidal_Zones',
  fileFormat: 'SHP'});
  
Export.image.toDrive({
  image: zones, 
  description:'Zones_AWEI_Linguado',
  folder: 'Intertidal_Zones', 
  fileNamePrefix:'Zones_AWEI_Linguado',
  region:geometry, 
  scale: 10,
  maxPixels: 100e9});
*/
//// Exportar vetores
var task = Export.table.toDrive({
  collection: vectors,
  description:'Intertidal_Zones_aweiValeMixed',
  folder: 'Intertidal_Zones',
  fileFormat: 'SHP'});
  
//// Exportar imagens
var ExportCol = function (col, folder, scale, tp, maxPixels, region) {
  scale = 10,
  maxPixels =  100e9;
  
  var nimg=col.size().getInfo();
  var colList = col.toList(nimg);
  var n = colList.size().getInfo();

  for (var i = 0 ; i < n; i++) {
    var img = ee.Image(colList.get(i));
    var id = img.get('system:id').getInfo() //.id().getInfo();
    region = region; // img.geometry().bounds().getInfo()["coordinates"];

      var imgtype = {"float":img.toFloat(), 
                 "byte":img.toByte(), 
                 "int":img.toInt(),
                 "double":img.toDouble()};

      Export.image.toDrive({
        image:imgtype[tp],
        description: id,
        folder: folder,
        fileNamePrefix: id + "_ndwi", 
        crs: 'EPSG:4326',
        region: region,
        scale: scale,
        maxPixels: maxPixels}
        );
        
      Export.image.toAsset({
        image:imgtype[tp],
        description: 'asset_'+id,
        assetId: 'projects/ee-index-images/assets/valeMixed/'+id,
        crs:'EPSG:4326',
        region: region,
        scale: scale,
        maxPixels: maxPixels}
        )}
};
var task2 = ExportCol(intertidal_zones, 'Intertidal_Zones', 10, 'float', 100e9, geometry);


Export.image.toDrive({
  image: NWI_STD, 
  description:'Zones_AWEI_ValeMixed',
  folder: 'Intertidal_Zones', 
  crs: 'EPSG:4326',
  fileNamePrefix:'awei_std_ValeMud',
  region:geometry, 
  scale: 10,
  maxPixels: 100e9});

///////////////////////// ==== testes ==== /////////////////////////

Map.addLayer(intertidal_zones.first(), vis_params, 'intertidal zone - first');
var imgList = intertidal_zones.toList(intertidal_zones.size());
Map.addLayer(ee.Image(imgList.get(-1)), vis_params, 'intertidal zone - last');
print("intertidal zones collection", intertidal_zones);