var startDate = ee.Date("2023-7-20");
var endDate = ee.Date("2023-8-20");

//bands of danyland
var dfBands = [
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
  "DOY",
  "stage",
];
//bands of paddy
var pfBands = [
  "B3",
  "B4",
  "B5",
  "B6",
  "B7",
  "B8",
  "B8A",
  "SZA",
  "DOY",
  "stage",
];

var laiVisualizationVis = {
  min: 0.0,
  max: 6.0,
  palette: ["e1e4b4", "999d60", "2ec409", "0a4b06"],
};

//crop mask
var cropMap = lulc.lt(3);
var dfMap = cropMap.gt(0);
var pfMap = lulc.lt(1);

var imRegion = region.filter(ee.Filter.eq("id", 5));
var spRegion = region.filter(ee.Filter.eq("id", 2));
var lkRegion = region.filter(ee.Filter.eq("id", 6));
var slRegion = region.filter(ee.Filter.eq("id", 4));
var lpRegion = region.filter(ee.Filter.eq("id", 1));

var s2 = ee.ImageCollection("COPERNICUS/S2_SR").filterDate(startDate, endDate);

var imImages = s2.filterBounds(imRegion);
var spImages = s2.filterBounds(spRegion);
var lkImages = s2.filterBounds(lkRegion);
var slImages = s2.filterBounds(slRegion);
var lpImages = s2.filterBounds(lpRegion);
var pfimages = s2.filterBounds(region);

//phenology feature
var imStages = [138, 191, 216, 249, 268, 311];
var spStages = [149, 183, 198, 243, 264, 307];
var lkStages = [139, 184, 202, 245, 264, 311];
var slStages = [144, 182, 199, 245, 263, 302];
var lpStages = [137, 184, 205, 241, 266, 319];
var pfStages = [149, 180, 193, 249, 271, 318];

var imLaiMap = imImages
  .map(preProcessImage)
  .map(getStageIM)
  .select(dfBands)
  .map(scaleDataIM)
  .map(calculateImLai)
  .max()
  .clip(imRegion)
  .updateMask(dfMap);

var spLaiMap = spImages
  .map(preProcessImage)
  .map(getStageSP)
  .select(dfBands)
  .map(scaleDataSP)
  .map(calculateSpLai)
  .max()
  .clip(spRegion)
  .updateMask(dfMap);

var lkLaiMap = lkImages
  .map(preProcessImage)
  .map(getStageLK)
  .select(dfBands)
  .map(scaleDataLK)
  .map(calculateLkLai)
  .max()
  .clip(lkRegion)
  .updateMask(dfMap);

var slLaiMap = slImages
  .map(preProcessImage)
  .map(getStageSL)
  .select(dfBands)
  .map(scaleDataSL)
  .map(calculateSlLai)
  .max()
  .clip(slRegion)
  .updateMask(dfMap);

var lpLaiMap = lpImages
  .map(preProcessImage)
  .map(getStageLP)
  .select(dfBands)
  .map(scaleDataLP)
  .map(calculateLpLai)
  .max()
  .clip(lpRegion)
  .updateMask(dfMap);

var pfLaiMap = pfimages
  .map(preProcessImage)
  .map(getStagePf)
  .select(pfBands)
  .map(scaleDataPf)
  .map(calculatePfLai)
  .max()
  .clip(region)
  .updateMask(pfMap);

var laiMaps = ee.ImageCollection([
  imLaiMap,
  spLaiMap,
  lkLaiMap,
  lpLaiMap,
  slLaiMap,
  pfLaiMap,
]);
var laiMap = laiMaps.mosaic();

Map.addLayer(laiMap, laiVisualizationVis, "DynoCrop");
Map.setCenter(131.75, 46.767, 11);

function preProcessImage(image) {
  var dayOfYear = image.date().getRelative("day", "year");
  var doyImage = ee.Image(dayOfYear).rename("DOY");
  var solarZenithAngle = getSolarZenithAngle(image);
  var scaledImage = scaleImage(image);
  var cloudlessImage = filtercloud(scaledImage);
  var snowlessImage = maskSnow(cloudlessImage);
  var result = snowlessImage.addBands(solarZenithAngle).addBands(doyImage);
  var mask = result.select("B8").gt(0);
  return result.updateMask(mask).set("date", dayOfYear);
}

function getSolarZenithAngle(image) {
  var solarZenithAngle = ee.Number(image.get("MEAN_SOLAR_ZENITH_ANGLE"));
  var projection = ee.Image(image.select("B3")).projection();
  var szaImage = ee
    .Image(solarZenithAngle)
    .setDefaultProjection(projection)
    .toFloat()
    .rename("SZA");
  return szaImage;
}

function maskSnow(image) {
  var b3 = image.select("B3");
  var b11 = image.select("B11");
  var ndsi = b3.subtract(b11).divide(b3.add(b11));
  var mask = ndsi.lt(0.4);
  return image.updateMask(mask);
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

function getStageIM(image) {
  var SMIN = imStages[0];
  var EMIN = imStages[5];
  var SOS = imStages[1];
  var SOP = imStages[2];
  var EOP = imStages[3];
  var EOS = imStages[4];
  var doy = ee.Number(image.get("date"));
  var stage = ee.Algorithms.If(
    doy.lt(SMIN),
    1,
    ee.Algorithms.If(
      doy.gte(EMIN),
      2,
      ee.Algorithms.If(
        doy.gte(EOS).and(doy.lt(EMIN)),
        3,
        ee.Algorithms.If(
          doy.gte(SMIN).and(doy.lt(SOS)),
          4,
          ee.Algorithms.If(
            doy.gte(SOS).and(doy.lt(SOP)),
            5,
            ee.Algorithms.If(
              doy.gte(EOP).and(doy.lt(EOS)),
              6,
              ee.Algorithms.If(doy.gte(SOP).and(doy.lt(EOP)), 7, null)
            )
          )
        )
      )
    )
  );
  var stageImage = ee.Image.constant(ee.Number(stage)).rename("stage");
  return image.addBands(stageImage);
}

function getStageSP(image) {
  var SMIN = spStages[0];
  var EMIN = spStages[5];
  var SOS = spStages[1];
  var SOP = spStages[2];
  var EOP = spStages[3];
  var EOS = spStages[4];
  var doy = ee.Number(image.get("date"));
  var stage = ee.Algorithms.If(
    doy.lt(SMIN),
    1,
    ee.Algorithms.If(
      doy.gte(EMIN),
      2,
      ee.Algorithms.If(
        doy.gte(EOS).and(doy.lt(EMIN)),
        3,
        ee.Algorithms.If(
          doy.gte(SMIN).and(doy.lt(SOS)),
          4,
          ee.Algorithms.If(
            doy.gte(SOS).and(doy.lt(SOP)),
            5,
            ee.Algorithms.If(
              doy.gte(EOP).and(doy.lt(EOS)),
              6,
              ee.Algorithms.If(doy.gte(SOP).and(doy.lt(EOP)), 7, null)
            )
          )
        )
      )
    )
  );
  var stageImage = ee.Image.constant(ee.Number(stage)).rename("stage");
  return image.addBands(stageImage);
}

function getStageLK(image) {
  var SMIN = lkStages[0];
  var EMIN = lkStages[5];
  var SOS = lkStages[1];
  var SOP = lkStages[2];
  var EOP = lkStages[3];
  var EOS = lkStages[4];
  var doy = ee.Number(image.get("date"));
  var stage = ee.Algorithms.If(
    doy.lt(SMIN),
    1,
    ee.Algorithms.If(
      doy.gte(EMIN),
      2,
      ee.Algorithms.If(
        doy.gte(EOS).and(doy.lt(EMIN)),
        3,
        ee.Algorithms.If(
          doy.gte(SMIN).and(doy.lt(SOS)),
          4,
          ee.Algorithms.If(
            doy.gte(SOS).and(doy.lt(SOP)),
            5,
            ee.Algorithms.If(
              doy.gte(EOP).and(doy.lt(EOS)),
              6,
              ee.Algorithms.If(doy.gte(SOP).and(doy.lt(EOP)), 7, null)
            )
          )
        )
      )
    )
  );
  var stageImage = ee.Image.constant(ee.Number(stage)).rename("stage");
  return image.addBands(stageImage);
}

function getStageLP(image) {
  var SMIN = lpStages[0];
  var EMIN = lpStages[5];
  var SOS = lpStages[1];
  var SOP = lpStages[2];
  var EOP = lpStages[3];
  var EOS = lpStages[4];
  var doy = ee.Number(image.get("date"));
  var stage = ee.Algorithms.If(
    doy.lt(SMIN),
    1,
    ee.Algorithms.If(
      doy.gte(EMIN),
      2,
      ee.Algorithms.If(
        doy.gte(EOS).and(doy.lt(EMIN)),
        3,
        ee.Algorithms.If(
          doy.gte(SMIN).and(doy.lt(SOS)),
          4,
          ee.Algorithms.If(
            doy.gte(SOS).and(doy.lt(SOP)),
            5,
            ee.Algorithms.If(
              doy.gte(EOP).and(doy.lt(EOS)),
              6,
              ee.Algorithms.If(doy.gte(SOP).and(doy.lt(EOP)), 7, null)
            )
          )
        )
      )
    )
  );
  var stageImage = ee.Image.constant(ee.Number(stage)).rename("stage");
  return image.addBands(stageImage);
}

function getStageSL(image) {
  var SMIN = slStages[0];
  var EMIN = slStages[5];
  var SOS = slStages[1];
  var SOP = slStages[2];
  var EOP = slStages[3];
  var EOS = slStages[4];
  var doy = ee.Number(image.get("date"));
  var stage = ee.Algorithms.If(
    doy.lt(SMIN),
    1,
    ee.Algorithms.If(
      doy.gte(EMIN),
      2,
      ee.Algorithms.If(
        doy.gte(EOS).and(doy.lt(EMIN)),
        3,
        ee.Algorithms.If(
          doy.gte(SMIN).and(doy.lt(SOS)),
          4,
          ee.Algorithms.If(
            doy.gte(SOS).and(doy.lt(SOP)),
            5,
            ee.Algorithms.If(
              doy.gte(EOP).and(doy.lt(EOS)),
              6,
              ee.Algorithms.If(doy.gte(SOP).and(doy.lt(EOP)), 7, null)
            )
          )
        )
      )
    )
  );
  var stageImage = ee.Image.constant(ee.Number(stage)).rename("stage");
  return image.addBands(stageImage);
}

function getStagePf(image) {
  var SMIN = pfStages[0];
  var EMIN = pfStages[5];
  var SOS = pfStages[1];
  var SOP = pfStages[2];
  var EOP = pfStages[3];
  var EOS = pfStages[4];
  var doy = ee.Number(image.get("date"));
  var stage = ee.Algorithms.If(
    doy.lt(SMIN),
    1,
    ee.Algorithms.If(
      doy.gte(EMIN),
      2,
      ee.Algorithms.If(
        doy.gte(EOS).and(doy.lt(EMIN)),
        3,
        ee.Algorithms.If(
          doy.gte(SMIN).and(doy.lt(SOS)),
          4,
          ee.Algorithms.If(
            doy.gte(SOS).and(doy.lt(SOP)),
            5,
            ee.Algorithms.If(
              doy.gte(EOP).and(doy.lt(EOS)),
              6,
              ee.Algorithms.If(doy.gte(SOP).and(doy.lt(EOP)), 7, null)
            )
          )
        )
      )
    )
  );
  var stageImage = ee.Image.constant(stage).rename("stage");
  return image.addBands(stageImage);
}

function scaleDataIM(image) {
  var list1 = [
    0.12921045, 0.14681306, 0.18009175, 0.2780241, 0.33931169, 0.35352299,
    0.36367761, 0.26667066, 0.19674214, 39.53298041, 222.90861045, 4.44836608,
  ];

  var list2 = [
    0.08141841, 0.10128002, 0.10025898, 0.06928274, 0.08447724, 0.07989593,
    0.08176071, 0.11111016, 0.11590532, 13.35135151, 69.39787412, 1.99849077,
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
  var SZA_mean = ee.Number(list1[9]);
  var DOY_mean = ee.Number(list1[10]);
  var stage_mean = ee.Number(list1[11]);

  var B3_std = ee.Number(list2[0]);
  var B4_std = ee.Number(list2[1]);
  var B5_std = ee.Number(list2[2]);
  var B6_std = ee.Number(list2[3]);
  var B7_std = ee.Number(list2[4]);
  var B8_std = ee.Number(list2[5]);
  var B8a_std = ee.Number(list2[6]);
  var B11_std = ee.Number(list2[7]);
  var B12_std = ee.Number(list2[8]);
  var SZA_std = ee.Number(list2[9]);
  var DOY_std = ee.Number(list2[10]);
  var stage_std = ee.Number(list2[11]);

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
  var SZA = image
    .select("SZA")
    .subtract(ee.Image(SZA_mean))
    .divide(ee.Image(SZA_std))
    .rename("SZA");
  var DOY = image
    .select("DOY")
    .subtract(ee.Image(DOY_mean))
    .divide(ee.Image(DOY_std))
    .rename("DOY");
  var stage = image
    .select("stage")
    .subtract(ee.Image(stage_mean))
    .divide(ee.Image(stage_std))
    .rename("stage");
  return B3.addBands(B4)
    .addBands(B5)
    .addBands(B6)
    .addBands(B7)
    .addBands(B8)
    .addBands(B8A)
    .addBands(B11)
    .addBands(B12)
    .addBands(SZA)
    .addBands(DOY)
    .addBands(stage);
}

function scaleDataSP(image) {
  var list1 = [
    0.09058886, 0.08552647, 0.12076861, 0.24925419, 0.31575583, 0.32503868,
    0.33643363, 0.20958293, 0.13704012, 39.21094161, 224.12166996, 5.12289372,
  ];

  var list2 = [
    0.06137901, 0.07143524, 0.0688106, 0.0802283, 0.11210017, 0.10735873,
    0.11173931, 0.06492216, 0.07893467, 10.9217346, 51.5355348, 2.01563741,
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
  var SZA_mean = ee.Number(list1[9]);
  var DOY_mean = ee.Number(list1[10]);
  var stage_mean = ee.Number(list1[11]);

  var B3_std = ee.Number(list2[0]);
  var B4_std = ee.Number(list2[1]);
  var B5_std = ee.Number(list2[2]);
  var B6_std = ee.Number(list2[3]);
  var B7_std = ee.Number(list2[4]);
  var B8_std = ee.Number(list2[5]);
  var B8a_std = ee.Number(list2[6]);
  var B11_std = ee.Number(list2[7]);
  var B12_std = ee.Number(list2[8]);
  var SZA_std = ee.Number(list2[9]);
  var DOY_std = ee.Number(list2[10]);
  var stage_std = ee.Number(list2[11]);

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
  var SZA = image
    .select("SZA")
    .subtract(ee.Image(SZA_mean))
    .divide(ee.Image(SZA_std))
    .rename("SZA");
  var DOY = image
    .select("DOY")
    .subtract(ee.Image(DOY_mean))
    .divide(ee.Image(DOY_std))
    .rename("DOY");
  var stage = image
    .select("stage")
    .subtract(ee.Image(stage_mean))
    .divide(ee.Image(stage_std))
    .rename("stage");
  return B3.addBands(B4)
    .addBands(B5)
    .addBands(B6)
    .addBands(B7)
    .addBands(B8)
    .addBands(B8A)
    .addBands(B11)
    .addBands(B12)
    .addBands(SZA)
    .addBands(DOY)
    .addBands(stage);
}

function scaleDataSL(image) {
  var list1 = [
    0.09606879, 0.10105468, 0.13314419, 0.24216094, 0.30608804, 0.31613503,
    0.3286448, 0.23471469, 0.16512905, 38.2239028, 221.07632707, 4.6497377,
  ];

  var list2 = [
    0.0454957, 0.06633366, 0.06466057, 0.0594656, 0.09416046, 0.08651523,
    0.09231606, 0.08702417, 0.10128751, 12.04761715, 57.42940868, 1.98664363,
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
  var SZA_mean = ee.Number(list1[9]);
  var DOY_mean = ee.Number(list1[10]);
  var stage_mean = ee.Number(list1[11]);

  var B3_std = ee.Number(list2[0]);
  var B4_std = ee.Number(list2[1]);
  var B5_std = ee.Number(list2[2]);
  var B6_std = ee.Number(list2[3]);
  var B7_std = ee.Number(list2[4]);
  var B8_std = ee.Number(list2[5]);
  var B8a_std = ee.Number(list2[6]);
  var B11_std = ee.Number(list2[7]);
  var B12_std = ee.Number(list2[8]);
  var SZA_std = ee.Number(list2[9]);
  var DOY_std = ee.Number(list2[10]);
  var stage_std = ee.Number(list2[11]);

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
  var SZA = image
    .select("SZA")
    .subtract(ee.Image(SZA_mean))
    .divide(ee.Image(SZA_std))
    .rename("SZA");
  var DOY = image
    .select("DOY")
    .subtract(ee.Image(DOY_mean))
    .divide(ee.Image(DOY_std))
    .rename("DOY");
  var stage = image
    .select("stage")
    .subtract(ee.Image(stage_mean))
    .divide(ee.Image(stage_std))
    .rename("stage");
  return B3.addBands(B4)
    .addBands(B5)
    .addBands(B6)
    .addBands(B7)
    .addBands(B8)
    .addBands(B8A)
    .addBands(B11)
    .addBands(B12)
    .addBands(SZA)
    .addBands(DOY)
    .addBands(stage);
}

function scaleDataLK(image) {
  var list1 = [
    0.10253944, 0.10950876, 0.14146267, 0.24092697, 0.29726513, 0.3083374,
    0.31912825, 0.23423075, 0.16685132, 39.56605495, 220.30859799, 4.50073858,
  ];

  var list2 = [
    0.05111164, 0.06975573, 0.06735431, 0.0681869, 0.10091497, 0.09458946,
    0.09999578, 0.08049443, 0.09074127, 12.94788819, 64.47867613, 2.07152838,
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
  var SZA_mean = ee.Number(list1[9]);
  var DOY_mean = ee.Number(list1[10]);
  var stage_mean = ee.Number(list1[11]);

  var B3_std = ee.Number(list2[0]);
  var B4_std = ee.Number(list2[1]);
  var B5_std = ee.Number(list2[2]);
  var B6_std = ee.Number(list2[3]);
  var B7_std = ee.Number(list2[4]);
  var B8_std = ee.Number(list2[5]);
  var B8a_std = ee.Number(list2[6]);
  var B11_std = ee.Number(list2[7]);
  var B12_std = ee.Number(list2[8]);
  var SZA_std = ee.Number(list2[9]);
  var DOY_std = ee.Number(list2[10]);
  var stage_std = ee.Number(list2[11]);

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
  var SZA = image
    .select("SZA")
    .subtract(ee.Image(SZA_mean))
    .divide(ee.Image(SZA_std))
    .rename("SZA");
  var DOY = image
    .select("DOY")
    .subtract(ee.Image(DOY_mean))
    .divide(ee.Image(DOY_std))
    .rename("DOY");
  var stage = image
    .select("stage")
    .subtract(ee.Image(stage_mean))
    .divide(ee.Image(stage_std))
    .rename("stage");
  return B3.addBands(B4)
    .addBands(B5)
    .addBands(B6)
    .addBands(B7)
    .addBands(B8)
    .addBands(B8A)
    .addBands(B11)
    .addBands(B12)
    .addBands(SZA)
    .addBands(DOY)
    .addBands(stage);
}

function scaleDataLP(image) {
  var list1 = [
    0.11753042, 0.12585969, 0.16016577, 0.26559281, 0.32886397, 0.33603034,
    0.34834243, 0.23867306, 0.17327024, 36.10195104, 223.09819921, 4.57934721,
  ];

  var list2 = [
    0.07253265, 0.08810723, 0.08599597, 0.0725575, 0.09806691, 0.09009409,
    0.09616253, 0.08570034, 0.10183844, 13.46483646, 65.92499859, 1.85429615,
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
  var SZA_mean = ee.Number(list1[9]);
  var DOY_mean = ee.Number(list1[10]);
  var stage_mean = ee.Number(list1[11]);

  var B3_std = ee.Number(list2[0]);
  var B4_std = ee.Number(list2[1]);
  var B5_std = ee.Number(list2[2]);
  var B6_std = ee.Number(list2[3]);
  var B7_std = ee.Number(list2[4]);
  var B8_std = ee.Number(list2[5]);
  var B8a_std = ee.Number(list2[6]);
  var B11_std = ee.Number(list2[7]);
  var B12_std = ee.Number(list2[8]);
  var SZA_std = ee.Number(list2[9]);
  var DOY_std = ee.Number(list2[10]);
  var stage_std = ee.Number(list2[11]);

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
  var SZA = image
    .select("SZA")
    .subtract(ee.Image(SZA_mean))
    .divide(ee.Image(SZA_std))
    .rename("SZA");
  var DOY = image
    .select("DOY")
    .subtract(ee.Image(DOY_mean))
    .divide(ee.Image(DOY_std))
    .rename("DOY");
  var stage = image
    .select("stage")
    .subtract(ee.Image(stage_mean))
    .divide(ee.Image(stage_std))
    .rename("stage");
  return B3.addBands(B4)
    .addBands(B5)
    .addBands(B6)
    .addBands(B7)
    .addBands(B8)
    .addBands(B8A)
    .addBands(B11)
    .addBands(B12)
    .addBands(SZA)
    .addBands(DOY)
    .addBands(stage);
}

function scaleDataPf(image) {
  var list1 = [
    0.10791229, 0.09912042, 0.16060573, 0.27486615, 0.32846314, 0.34057018,
    0.35264059, 39.60527981, 235.66208302, 5.55752567,
  ];

  var list2 = [
    0.05518398, 0.06373726, 0.06092403, 0.06305802, 0.08083667, 0.0795217,
    0.08368158, 10.1219361, 46.31829794, 1.6380607,
  ];

  var B3_mean = ee.Number(list1[0]);
  var B4_mean = ee.Number(list1[1]);
  var B5_mean = ee.Number(list1[2]);
  var B6_mean = ee.Number(list1[3]);
  var B7_mean = ee.Number(list1[4]);
  var B8_mean = ee.Number(list1[5]);
  var B8a_mean = ee.Number(list1[6]);
  var SZA_mean = ee.Number(list1[7]);
  var DOY_mean = ee.Number(list1[8]);
  var stage_mean = ee.Number(list1[9]);

  var B3_std = ee.Number(list2[0]);
  var B4_std = ee.Number(list2[1]);
  var B5_std = ee.Number(list2[2]);
  var B6_std = ee.Number(list2[3]);
  var B7_std = ee.Number(list2[4]);
  var B8_std = ee.Number(list2[5]);
  var B8a_std = ee.Number(list2[6]);
  var SZA_std = ee.Number(list2[7]);
  var DOY_std = ee.Number(list2[8]);
  var stage_std = ee.Number(list2[9]);

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
  var SZA = image
    .select("SZA")
    .subtract(ee.Image(SZA_mean))
    .divide(ee.Image(SZA_std))
    .rename("SZA");
  var DOY = image
    .select("DOY")
    .subtract(ee.Image(DOY_mean))
    .divide(ee.Image(DOY_std))
    .rename("DOY");
  var stage = image
    .select("stage")
    .subtract(ee.Image(stage_mean))
    .divide(ee.Image(stage_std))
    .rename("stage");
  return B3.addBands(B4)
    .addBands(B5)
    .addBands(B6)
    .addBands(B7)
    .addBands(B8)
    .addBands(B8A)
    .addBands(SZA)
    .addBands(DOY)
    .addBands(stage);
}

function calculateSpLai(image) {
  var w1 = ee.Array([
    [0.34762806, -0.15203889, 0.21508104, -0.39013442],
    [0.15197675, -0.07951006, 0.52898586, 0.02490202],
    [0.17642346, -0.21441214, 0.3251366, 0.13444504],
    [-0.2622631, 0.2192583, 0.09400921, -0.03752564],
    [0.16336621, 0.20212737, -0.49275428, -0.738163],
    [-0.50123465, 0.571135, -0.4690471, -0.400031],
    [-0.39362082, 0.47111684, -0.21595481, 0.25611782],
    [0.12192989, -0.83451587, -0.22162186, 0.18262246],
    [-0.093978, -0.78531486, 0.22714022, -0.00318635],
    [-0.06342559, -0.80566484, 0.04543541, -0.16140383],
    [-0.4206227, 0.48375213, -0.36887878, 0.0629869],
    [0.57843256, 0.55556184, 0.02084534, 0.12010793],
  ]);
  var b1 = ee.Array([[0.89449435, 0.5092986, -0.45082048, -0.2627227]]);
  var w2 = ee.Array([[0.57363564], [0.96167916], [-0.63058496], [-0.26699096]]);
  var b2 = ee.Number(0.4727782);
  var w1Image = ee.Image(w1).toArray().toArray(1);
  var w2Image = ee.Image(w2).toArray().toArray(1);
  var b1Image = ee.Image(b1).toArray();
  var imageArray = image.toArray().toArray(1).matrixTranspose();
  var w1ImageArray = imageArray.matrixMultiply(w1Image).add(b1Image);
  var w2ImageArray = w1ImageArray.abs().add(w1ImageArray).divide(2);
  var lai = w2ImageArray.matrixMultiply(w2Image).add(b2).matrixDeterminant();
  return ee.Image(lai);
}

function calculateImLai(image) {
  var w1 = ee.Array([
    [0.24858776, 0.6112801, 0.12482149, 0.04555824],
    [0.39682665, -0.4356552, -0.338933, -0.00303169],
    [0.30234838, -0.29420143, -0.5263643, 0.21601048],
    [-0.00081019, 0.35535866, 0.34451446, 0.23722158],
    [0.2640344, 0.23429346, 0.4343099, -0.5120991],
    [-0.03896748, 0.60618526, -0.2769492, 0.25939977],
    [0.18943913, -0.06687357, -0.02566948, -0.3680175],
    [-0.21798068, -0.85252184, -0.8706433, -0.32852736],
    [-0.41342807, -0.28949258, -0.10179424, 0.39927477],
    [0.13409947, -1.1278452, -0.2236874, 0.01474985],
    [0.35347855, -0.7112764, -0.11455248, -0.06974268],
    [-0.8394318, 0.7940371, -0.08989486, -0.531798],
  ]);
  var b1 = ee.Array([[-0.26651904, 0.17211622, 0.45587882, -0.44848228]]);
  var w2 = ee.Array([[-0.17579319], [0.35917473], [0.5214407], [-0.3303467]]);
  var b2 = ee.Number(0.4631008);
  var w1Image = ee.Image(w1).toArray().toArray(1);
  var w2Image = ee.Image(w2).toArray().toArray(1);
  var b1Image = ee.Image(b1).toArray();
  var imageArray = image.toArray().toArray(1).matrixTranspose();
  var w1ImageArray = imageArray.matrixMultiply(w1Image).add(b1Image);
  var w2ImageArray = w1ImageArray.abs().add(w1ImageArray).divide(2);
  var lai = w2ImageArray.matrixMultiply(w2Image).add(b2).matrixDeterminant();
  return ee.Image(lai);
}

function calculateSlLai(image) {
  var w1 = ee.Array([
    [0.56139505, 0.19134074, -0.5939529, 0.41583043],
    [-0.3593119, 0.5144039, -0.19915232, -0.2017282],
    [-0.26872256, 0.07616501, -0.73628426, 0.03257864],
    [0.5167542, 0.43944827, -0.46775806, 0.61196816],
    [-0.23533265, -0.28673142, 0.3326161, -0.393642],
    [-0.0287088, 0.20681801, 0.49581736, 0.03377418],
    [-0.30503273, -0.5345772, 0.3613783, -0.20347239],
    [0.06377839, 0.29985526, 0.14908247, 0.2564752],
    [-0.5222149, 0.21704072, -0.11927287, 0.20433849],
    [-0.41213366, -0.3932473, 0.0680889, -0.2725923],
    [0.0140965, -0.23066388, -0.09708909, -0.20643526],
    [0.8331336, -0.29439467, -0.04700631, 0.12666607],
  ]);
  var b1 = ee.Array([[0.33513376, -0.16382405, 0.2218325, 0.44143823]]);
  var w2 = ee.Array([[0.9115111], [-0.47005185], [0.7863187], [0.6128489]]);
  var b2 = ee.Number(0.306371);
  var w1Image = ee.Image(w1).toArray().toArray(1);
  var w2Image = ee.Image(w2).toArray().toArray(1);
  var b1Image = ee.Image(b1).toArray();
  var imageArray = image.toArray().toArray(1).matrixTranspose();
  var w1ImageArray = imageArray.matrixMultiply(w1Image).add(b1Image);
  var w2ImageArray = w1ImageArray.abs().add(w1ImageArray).divide(2);
  var lai = w2ImageArray.matrixMultiply(w2Image).add(b2).matrixDeterminant();
  return ee.Image(lai);
}

function calculateLpLai(image) {
  var w1 = ee.Array([
    [0.6425398, -0.26521665, -0.5724022, 0.30541033],
    [-0.20696153, 0.2670934, 0.33000782, -0.18930234],
    [-0.61326003, -0.11386528, 0.39018995, 0.00594914],
    [0.6099209, -0.01745109, -0.08025566, -0.42255563],
    [0.44455722, 0.09729426, -0.09467651, -0.2107149],
    [-0.13463528, -0.29746693, 0.29275572, -0.10236961],
    [0.41581193, 0.16855694, -0.41056693, -0.15456317],
    [-0.17616002, -0.4400503, 0.24119107, 0.33992177],
    [-0.41220585, 0.11498795, -0.45544901, -0.163587],
    [0.26946846, -0.6222021, -0.5932897, -0.4593545],
    [0.07648234, -0.20547187, 0.30484378, 0.03941961],
    [0.04463359, 0.7196938, 0.200971, 0.4901846],
  ]);
  var b1 = ee.Array([[0.35734388, 0.23873986, -0.17643623, -0.2021954]]);
  var w2 = ee.Array([[0.3079062], [1.2474533], [-0.45446223], [-0.6385144]]);
  var b2 = ee.Number(0.2005389);
  var w1Image = ee.Image(w1).toArray().toArray(1);
  var w2Image = ee.Image(w2).toArray().toArray(1);
  var b1Image = ee.Image(b1).toArray();
  var imageArray = image.toArray().toArray(1).matrixTranspose();
  var w1ImageArray = imageArray.matrixMultiply(w1Image).add(b1Image);
  var w2ImageArray = w1ImageArray.abs().add(w1ImageArray).divide(2);
  var lai = w2ImageArray.matrixMultiply(w2Image).add(b2).matrixDeterminant();
  return ee.Image(lai);
}

function calculateLkLai(image) {
  var w1 = ee.Array([
    [-0.00664338, -0.1864975, -0.6019374, -0.28107795],
    [-0.14512597, -0.4672227, 0.0230937, 0.42901832],
    [0.27221173, -0.4960425, 0.37939465, 0.37523708],
    [0.00100587, 0.08813483, 0.06411487, -0.22918189],
    [0.664333, 0.76893276, 0.4271323, 0.10383198],
    [-0.29164103, -0.3152204, -0.02597727, 0.50172305],
    [-0.40261102, 0.22685918, 0.00423189, 0.08260678],
    [-0.23212913, 0.16729859, -0.07069009, -0.48205695],
    [0.30717367, 0.21709424, -0.12732811, -0.7870764],
    [0.10635711, -0.49172527, -0.0465995, -0.17719896],
    [0.00709331, 0.09497648, 0.3244545, -0.31034067],
    [-0.39602357, -0.04762499, -0.0997052, 0.81465137],
  ]);
  var b1 = ee.Array([[-0.26935536, 0.31426826, -0.16965422, 0.36770687]]);
  var w2 = ee.Array([[-0.5597108], [0.9647082], [-0.36611083], [0.86536795]]);
  var b2 = ee.Number(0.34272093);
  var w1Image = ee.Image(w1).toArray().toArray(1);
  var w2Image = ee.Image(w2).toArray().toArray(1);
  var b1Image = ee.Image(b1).toArray();
  var imageArray = image.toArray().toArray(1).matrixTranspose();
  var w1ImageArray = imageArray.matrixMultiply(w1Image).add(b1Image);
  var w2ImageArray = w1ImageArray.abs().add(w1ImageArray).divide(2);
  var lai = w2ImageArray.matrixMultiply(w2Image).add(b2).matrixDeterminant();
  return ee.Image(lai);
}

function calculatePfLai(image) {
  var w1 = ee.Array([
    [-0.08580608, 0.05099415, -0.49858248, -0.09345081],
    [-0.20987754, 0.11611164, 0.6311445, 0.14016284],
    [0.58126104, 0.48490855, 0.43446603, -0.6520876],
    [0.15613021, -0.53248084, -0.66270185, 0.46619254],
    [-0.55537605, 0.00942937, 0.17108467, 0.37063447],
    [0.29502946, 0.23700312, -0.19784932, -0.38175637],
    [-0.32978162, -0.46679172, 0.23000735, 0.13513955],
    [0.3221178, 0.11806282, 0.3483637, -0.59266275],
    [0.67716336, 0.36694893, 0.9494433, -0.30303589],
    [-0.20658305, -0.36813655, -0.25891066, 0.7293286],
  ]);
  var b1 = ee.Array([[0.457575, -0.24695854, -0.365228, 0.39853358]]);
  var w2 = ee.Array([[0.74183536], [-0.24173722], [-0.55589545], [1.2039874]]);
  var b2 = ee.Number(0.43366498);
  var w1Image = ee.Image(w1).toArray().toArray(1);
  var w2Image = ee.Image(w2).toArray().toArray(1);
  var b1Image = ee.Image(b1).toArray();
  var imageArray = image.toArray().toArray(1).matrixTranspose();
  var w1ImageArray = imageArray.matrixMultiply(w1Image).add(b1Image);
  var w2ImageArray = w1ImageArray.abs().add(w1ImageArray).divide(2);
  var lai = w2ImageArray.matrixMultiply(w2Image).add(b2).matrixDeterminant();
  return ee.Image(lai);
}
