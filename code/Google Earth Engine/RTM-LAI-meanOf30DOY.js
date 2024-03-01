var startDate = ee.Date("2023-07-20");
var endDate = ee.Date("2023-08-20");

var s2 = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED");
var collection = s2.filterDate(startDate, endDate).filterBounds(region);

var proj = ee.Image(collection.first().select("B8")).projection();

var lulc = lulc.lt(3);

var lais = collection
  .map(preProcessImage)
  .map(scaleData)
  .map(calculateLai)
  .mean()
  .updateMask(lulc);

var laiVisualizationVis = {
  min: 0.0,
  max: 6.0,
  palette: ["e1e4b4", "999d60", "2ec409", "0a4b06"],
};

Map.addLayer(lais, laiVisualizationVis, "lai");
Map.setCenter(131.75, 46.767, 11);

function preProcessImage(image) {
  var date = image.date();
  var doy = date.getRelative("day", "year");
  var geometryInfo = getGeometryInfo(image);
  var scaledImage = scaleImage(image);
  var cloudlessImage = filtercloud(scaledImage);
  var cloudsnowlessImage = maskSnow(cloudlessImage);
  var result = cloudsnowlessImage.addBands(geometryInfo);
  var mask = result.select("B8").gt(0);
  return result.updateMask(mask).set("date", doy);
}

function getGeometryInfo(image) {
  var REA;
  var proj = ee.Image(image.select("B8")).projection();
  var SZA = ee.Number(image.get("MEAN_SOLAR_ZENITH_ANGLE"));
  var VZA = ee.Number(image.get("MEAN_INCIDENCE_ZENITH_ANGLE_B8"));
  var SAA = ee.Number(image.get("MEAN_SOLAR_AZIMUTH_ANGLE"));
  var VAA = ee.Number(image.get("MEAN_INCIDENCE_AZIMUTH_ANGLE_B8"));
  var rea = SAA.subtract(VAA).abs();
  REA = ee.Number(
    ee.Algorithms.If(rea.lt(ee.Number(180)), rea, ee.Number(360).subtract(rea))
  );
  var img_SZA = ee.Image.constant(SZA)
    .setDefaultProjection(proj)
    .toFloat()
    .rename("SZA");
  var img_VZA = ee.Image.constant(VZA)
    .setDefaultProjection(proj)
    .toFloat()
    .rename("VZA");
  var img_REA = ee.Image.constant(REA)
    .setDefaultProjection(proj)
    .toFloat()
    .rename("REA");
  var result = img_SZA.addBands(img_VZA).addBands(img_REA);
  return result;
}

function cloudScore(img) {
  var rescale = function (img, exp, thresholds) {
    return img
      .expression(exp, { img: img })
      .subtract(thresholds[0])
      .divide(thresholds[1] - thresholds[0]);
  };
  var score = ee.Image.constant(1.0);
  score = score.min(rescale(img, "img.blue", [0.1, 0.3]));
  score = score.min(rescale(img, "img.red + img.green + img.blue", [0.2, 0.8]));
  score = score.min(
    rescale(img, "img.nir + img.swir1 + img.swir2", [0.3, 0.8])
  );
  var ndsi = img.normalizedDifference(["green", "swir1"]);
  return score.min(rescale(ndsi, "img", [0.8, 0.6]));
}

function rmCloudByScore(image, thread) {
  var preBands = ["B2", "B3", "B4", "B8", "B11", "B12"];
  var newBands = ["blue", "green", "red", "nir", "swir1", "swir2"];
  var score = cloudScore(image.select(preBands, newBands));
  score = score.multiply(100).byte().rename("cloud");
  return image.addBands(score).updateMask(score.lte(thread));
}

function filtercloud(image) {
  return rmCloudByScore(image, 30);
}

function maskSnow(image) {
  var b3 = image.select("B3");
  var b11 = image.select("B11");
  var ndsi = b3.subtract(b11).divide(b3.add(b11));
  var mask = ndsi.lt(0.4);
  return image.updateMask(mask);
}

function scaleImage(image) {
  return ee
    .Image(
      image.select([
        "B2",
        "B3",
        "B4",
        "B5",
        "B6",
        "B7",
        "B8",
        "B8A",
        "B11",
        "B12",
      ])
    )
    .divide(10000);
}

function scaleData(image) {
  var list1 = [
    0.05378174, 0.03251672, 0.08764893, 0.31300118, 0.38789542, 0.39126857,
    0.39255144, 0.17181878, 0.07619735, 27.48212026, 5.98005043, 90.1514952,
  ];

  var list2 = [
    0.02425916, 0.02946091, 0.03129181, 0.07178698, 0.10299614, 0.10003491,
    0.09789412, 0.05614473, 0.05511926, 4.31597605, 2.31520729, 34.74916577,
  ];
  var B3_mean = ee.Number(list1[0]);
  var B4_mean = ee.Number(list1[1]);
  var B5_mean = ee.Number(list1[2]);
  var B6_mean = ee.Number(list1[3]);
  var B7_mean = ee.Number(list1[4]);
  var B8_mean = ee.Number(list1[5]);
  var B8a_mean = ee.Number(list1[6]);
  var B11_mean = ee.Number(list1[7]);
  var B12_mean = ee.Number(list1[8]);
  var sza_mean = ee.Number(list1[9]);
  var vza_mean = ee.Number(list1[10]);
  var rea_mean = ee.Number(list1[11]);

  var B3_std = ee.Number(list2[0]);
  var B4_std = ee.Number(list2[1]);
  var B5_std = ee.Number(list2[2]);
  var B6_std = ee.Number(list2[3]);
  var B7_std = ee.Number(list2[4]);
  var B8_std = ee.Number(list2[5]);
  var B8a_std = ee.Number(list2[6]);
  var B11_std = ee.Number(list2[7]);
  var B12_std = ee.Number(list2[8]);
  var sza_std = ee.Number(list2[9]);
  var vza_std = ee.Number(list2[10]);
  var rea_std = ee.Number(list2[11]);

  var B3 = image
    .select("B3")
    .subtract(ee.Image(B3_mean))
    .divide(ee.Image(B3_std))
    .rename("B3");
  var B4 = image
    .select("B4")
    .subtract(ee.Image(B4_mean))
    .divide(ee.Image(B4_std))
    .rename("B4");
  var B5 = image
    .select("B5")
    .subtract(ee.Image(B5_mean))
    .divide(ee.Image(B5_std))
    .rename("B5");
  var B6 = image
    .select("B6")
    .subtract(ee.Image(B6_mean))
    .divide(ee.Image(B6_std))
    .rename("B6");
  var B7 = image
    .select("B7")
    .subtract(ee.Image(B7_mean))
    .divide(ee.Image(B7_std))
    .rename("B7");
  var B8 = image
    .select("B8")
    .subtract(ee.Image(B8_mean))
    .divide(ee.Image(B8_std))
    .rename("B8");
  var B8A = image
    .select("B8A")
    .subtract(ee.Image(B8a_mean))
    .divide(ee.Image(B8a_std))
    .rename("B8A");
  var B11 = image
    .select("B11")
    .subtract(ee.Image(B11_mean))
    .divide(ee.Image(B11_std))
    .rename("B11");
  var B12 = image
    .select("B12")
    .subtract(ee.Image(B12_mean))
    .divide(ee.Image(B12_std))
    .rename("B12");
  var sza = image
    .select("SZA")
    .subtract(ee.Image(sza_mean))
    .divide(ee.Image(sza_std))
    .rename("SZA");
  var vza = image
    .select("VZA")
    .subtract(ee.Image(vza_mean))
    .divide(ee.Image(vza_std))
    .rename("VZA");
  var rea = image
    .select("REA")
    .subtract(ee.Image(rea_mean))
    .divide(ee.Image(rea_std))
    .rename("REA");
  var doy = image.get("date");
  return B3.addBands(B4)
    .addBands(B5)
    .addBands(B6)
    .addBands(B7)
    .addBands(B8)
    .addBands(B8A)
    .addBands(B11)
    .addBands(B12)
    .addBands(sza)
    .addBands(vza)
    .addBands(rea)
    .set("date", doy);
}

function calculateLai(image) {
  var e = ee.Number(Math.E);
  var eImage = ee.Image(e);

  var w1Layer = ee.Array([
    [-0.30189362, 0.04263008, -0.12200212, -0.5074681],
    [-0.31786388, -0.5385011, 0.11258744, 0.02595216],
    [-0.3065987, -0.14878222, 0.06267114, -0.74538565],
    [0.83257097, 0.8156011, 0.66078836, 0.21162407],
    [0.05875765, 0.19894247, 0.02305646, 0.8461102],
    [-0.15184069, 0.06088702, 0.5783646, -0.26660714],
    [0.4810504, 0.13298692, 0.3921261, 0.5832697],
    [-0.11354139, -0.6869614, -0.2896235, -0.2028416],
    [-0.1195024, -0.13483779, -0.788628, -0.4612247],
    [0.3704409, -0.27364618, -0.06532845, 0.03748444],
    [-0.5125734, 0.19796574, 0.03980537, 0.19192182],
    [0.03037837, -0.06963252, 0.20725842, -0.10221627],
  ]);

  var b1Layer = ee.Array([[0.32772833, 0.31321853, 0.2934478, 0.26315033]]);
  var w2Layer = ee.Array([[1.1414809], [1.4887092], [1.3530422], [1.4129509]]);
  var b2Layer = ee.Number(0.32621416);
  var w1LayerImage = ee.Image(w1Layer).toArray().toArray(1);
  var w2LayerImage = ee.Image(w2Layer).toArray().toArray(1);
  var b1LayerImage = ee.Image(b1Layer).toArray();
  var s2 = image.toArray().toArray(1).matrixTranspose();
  var layer1 = s2.matrixMultiply(w1LayerImage).add(b1LayerImage);
  var layer2 = ee
    .Image(1)
    .divide(eImage.pow(layer1.multiply(-1)).add(ee.Number(1)));
  var lai = layer2
    .matrixMultiply(w2LayerImage)
    .add(b2Layer)
    .matrixDeterminant();
  return lai;
}
