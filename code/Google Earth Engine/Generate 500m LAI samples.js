var startDate = ee.Date("2019-1-01");
var endDate = ee.Date("2019-4-20");
var lulc = lulc2019;
var id = 1;
var cloud = 50;
var roi = region.filter(ee.Filter.eq("id", id));

var bands = [
  "latitude",
  "longitude",
  "lai",
  "label",
  "B3",
  "B4",
  "B5",
  "B6",
  "B7",
  "B8",
  "B8A",
  "B11",
  "B12",
  "SZA",
  "VZA",
  "REA",
  "DOY",
];

var lai = ee.ImageCollection("MODIS/006/MCD15A3H").filterBounds(roi);
var s2 = ee.ImageCollection("COPERNICUS/S2_SR");

var collection_1 = s2
  .filterDate(startDate, endDate)
  .filterBounds(roi)
  .filter(ee.Filter.lt("CLOUD_COVERAGE_ASSESSMENT", cloud));

print(collection_1.size());

var s2Projection = ee.Image(collection_1.first().select("B8")).projection();
var s2Projection_500 = s2Projection.atScale(500);

var maxDiffFilter = ee.Filter.maxDifference({
  difference: 2.5 * 24 * 60 * 60 * 1000,
  leftField: "system:time_start",
  rightField: "system:time_start",
});

var saveBestJoin = ee.Join.saveBest({
  matchKey: "bestImage",
  measureKey: "timeDiff",
});

var overall = saveBestJoin.apply(collection_1, lai, maxDiffFilter);

var overall_list = overall.toList(collection_1.size());

function Laifpar(img) {
  var s2Projection = ee.Image(img.select("B8")).projection();
  var s2Projection_500 = s2Projection.atScale(500);
  var MCD15A3H = ee.Image(img.get("bestImage"));
  var lai = MCD15A3H.select("Lai");
  var mask_unsaturation = MCD15A3H.select("FparLai_QC")
    .bitwiseAnd(ee.Number.parse("11100000", 2))
    .eq(0);
  var mask_saturation = MCD15A3H.select("FparLai_QC")
    .bitwiseAnd(ee.Number.parse("11100000", 2))
    .eq(32);
  var unsaturation_label = mask_unsaturation.multiply(1);
  var saturation_label = mask_saturation.multiply(5);
  var label = unsaturation_label.add(saturation_label).rename("label");
  var lai_mask = label.gt(0);
  var lai_res = lai.updateMask(lai_mask).rename("lai");
  var res = lai_res.addBands(label);
  var lai_500 = res.reproject({ crs: s2Projection_500, scale: 500 });
  return lai_500;
}

function geometry(img) {
  var REA;
  var s2Projection = ee.Image(img.select("B8")).projection();
  var s2Projection_500 = s2Projection.atScale(500);
  var SZA = ee.Number(img.get("MEAN_SOLAR_ZENITH_ANGLE"));
  var VZAB2 = ee.Number(img.get("MEAN_INCIDENCE_ZENITH_ANGLE_B2"));
  var VZAB3 = ee.Number(img.get("MEAN_INCIDENCE_ZENITH_ANGLE_B3"));
  var VZAB4 = ee.Number(img.get("MEAN_INCIDENCE_ZENITH_ANGLE_B4"));
  var VZAB5 = ee.Number(img.get("MEAN_INCIDENCE_ZENITH_ANGLE_B5"));
  var VZAB6 = ee.Number(img.get("MEAN_INCIDENCE_ZENITH_ANGLE_B6"));
  var VZAB7 = ee.Number(img.get("MEAN_INCIDENCE_ZENITH_ANGLE_B7"));
  var VZAB8 = ee.Number(img.get("MEAN_INCIDENCE_ZENITH_ANGLE_B8"));
  var VZAB8A = ee.Number(img.get("MEAN_INCIDENCE_ZENITH_ANGLE_B8A"));
  var VZAB11 = ee.Number(img.get("MEAN_INCIDENCE_ZENITH_ANGLE_B11"));
  var VZAB12 = ee.Number(img.get("MEAN_INCIDENCE_ZENITH_ANGLE_B12"));
  var VZA = VZAB2.add(VZAB3)
    .add(VZAB4)
    .add(VZAB5)
    .add(VZAB6)
    .add(VZAB7)
    .add(VZAB8)
    .add(VZAB8A)
    .add(VZAB11)
    .add(VZAB12)
    .divide(10);
  var SAA = ee.Number(img.get("MEAN_SOLAR_AZIMUTH_ANGLE"));
  var VAAB2 = ee.Number(img.get("MEAN_INCIDENCE_AZIMUTH_ANGLE_B2"));
  var VAAB3 = ee.Number(img.get("MEAN_INCIDENCE_AZIMUTH_ANGLE_B3"));
  var VAAB4 = ee.Number(img.get("MEAN_INCIDENCE_AZIMUTH_ANGLE_B4"));
  var VAAB5 = ee.Number(img.get("MEAN_INCIDENCE_AZIMUTH_ANGLE_B5"));
  var VAAB6 = ee.Number(img.get("MEAN_INCIDENCE_AZIMUTH_ANGLE_B6"));
  var VAAB7 = ee.Number(img.get("MEAN_INCIDENCE_AZIMUTH_ANGLE_B7"));
  var VAAB8 = ee.Number(img.get("MEAN_INCIDENCE_AZIMUTH_ANGLE_B8"));
  var VAAB8A = ee.Number(img.get("MEAN_INCIDENCE_AZIMUTH_ANGLE_B8A"));
  var VAAB11 = ee.Number(img.get("MEAN_INCIDENCE_AZIMUTH_ANGLE_B11"));
  var VAAB12 = ee.Number(img.get("MEAN_INCIDENCE_AZIMUTH_ANGLE_B12"));
  var VAA = VAAB2.add(VAAB3)
    .add(VAAB4)
    .add(VAAB5)
    .add(VAAB6)
    .add(VAAB7)
    .add(VAAB8)
    .add(VAAB8A)
    .add(VAAB11)
    .add(VAAB12)
    .divide(10);
  var rea = SAA.subtract(VAA).abs();
  REA = ee.Number(
    ee.Algorithms.If(rea.lt(ee.Number(180)), rea, ee.Number(360).subtract(rea))
  );
  var img_SZA = ee
    .Image(SZA)
    .setDefaultProjection(s2Projection_500)
    .toFloat()
    .rename("SZA");
  var img_VZA = ee
    .Image(VZA)
    .setDefaultProjection(s2Projection_500)
    .toFloat()
    .rename("VZA");
  var img_REA = ee
    .Image(REA)
    .setDefaultProjection(s2Projection_500)
    .toFloat()
    .rename("REA");
  var res = img_SZA.addBands(img_VZA).addBands(img_REA);
  return img.addBands(res);
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

function resampleAndcvMask(img) {
  var s2Projection = ee.Image(img.select("B8")).projection();
  var s2Projection_500 = s2Projection.atScale(500);
  var mean = img
    .reduceResolution({
      reducer: ee.Reducer.mean(),
      maxPixels: 60000,
    })
    .reproject({ crs: s2Projection_500 });

  var mean_6 = mean.select(["B3", "B4", "B8"]);

  var std = img
    .reduceResolution({
      reducer: ee.Reducer.stdDev(),
      maxPixels: 60000,
    })
    .reproject({ crs: s2Projection_500 });
  var std_6 = std.select(["B3", "B4", "B8"]);
  var cv = std_6.divide(mean_6).rename(["cv3", "cv4", "cv8"]);
  var mask = cv
    .select("cv3")
    .lt(0.1)
    .multiply(cv.select("cv4").lt(0.1))
    .multiply(cv.select("cv8").lt(0.1));
  return mean.updateMask(mask);
}

//去云
var _cloudScore = function (img) {
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
};

function rmCloudByScore(image, thread) {
  var preBands = ["B2", "B3", "B4", "B8", "B11", "B12"];
  var newBands = ["blue", "green", "red", "nir", "swir1", "swir2"];
  var score = _cloudScore(image.select(preBands, newBands));
  score = score.multiply(100).byte().rename("cloud");
  return image.addBands(score).updateMask(score.lte(thread));
}

function Filtercloud(image) {
  return rmCloudByScore(image, 30);
}

function crop_mask(img) {
  var extent = img.geometry();
  var s2Projection = ee.Image(img.select("B8")).projection();
  var s2Projection_500 = s2Projection.atScale(500);
  var crop_30 = lulc.lt(3).clip(extent);
  var crop_500 = crop_30
    .reduceResolution({
      reducer: ee.Reducer.mean(),
      maxPixels: 60000,
    })
    .reproject({ crs: s2Projection_500 });
  return crop_500.gt(0.9);
}

function resample500(image) {
  var image_fiveBands = image.select(["B2", "B3", "B4", "B8", "B11", "B12"]);
  var image_fiveBands_500 = image_fiveBands
    .reduceResolution({
      reducer: ee.Reducer.mean(),
      maxPixels: 2600,
    })
    .reproject({ crs: s2Projection_500 });
  return image_fiveBands_500.divide(10000);
}

function preProcessing(img) {
  var s2Projection = ee.Image(img).select("B8").projection();
  var s2Projection_500 = s2Projection.atScale(500);
  img = ee.Image(img);
  var extent = img.geometry();
  var center = extent.centroid();
  var geo = geometry(img);
  var laifpar = Laifpar(img);
  var crop_500_mask = crop_mask(img);
  var Doy = ee
    .Image(ee.Number(img.date().getRelative("day", "year")))
    .rename("DOY");
  img = scaleImage(img);
  img = resampleAndcvMask(img);
  img = Filtercloud(img)
    .addBands(geo)
    .addBands(ee.Image.pixelLonLat())
    .addBands(Doy);
  img = img.updateMask(crop_500_mask);
  var selfmask = img.select("B8").gt(0).multiply(laifpar.select("label").gt(0));
  var res = img
    .addBands(laifpar)
    .updateMask(selfmask)
    .select(bands)
    .multiply(1000000)
    .toInt32()
    .clip(roi);

  var point = ee.Image(res).sample({
    region: extent,
    geometries: true,
    projection: s2Projection_500,
    scale: 500,
  });
  return point;
}

function exportToCsv(start, end, img_collection, bands, date) {
  Export.table.toDrive({
    description: start + "-" + end + date,
    collection: img_collection,
    selectors: bands,
  });
}

var collection_size = collection_1.size().toInt().getInfo();
var collection_size_s = collection_size / 10;

for (var i = 0; i < collection_size_s; i++) {
  var start = i * 10;
  var end = (i + 1) * 10;
  var overall_res = ee.ImageCollection(overall_list.slice(start, end, 1));
  var feature_collection = overall_res.map(preProcessing).flatten();
  print(start, end, feature_collection.size());
  exportToCsv(start, end, feature_collection, bands, "101-420");
}
