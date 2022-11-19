#ifndef _BIOMES_H_
#define _BIOMES_H_

#include <string>
#include <tuple>
#include "./vectorMath.h"

class UV {
public:
  float u;
  float v;
};

enum class BIOME : uint8_t {
  biOcean,
  biPlains,
  biDesert,
  biExtremeHills,
  biForest,
  biTaiga,
  biSwampland,
  biRiver,
  biNether,
  biEnd,
  biFrozenOcean,
  biFrozenRiver,
  biTundra,
  biIceMountains,
  biMushroomIsland,
  biMushroomShore,
  biBeach,
  biDesertHills,
  biForestHills,
  biTaigaHills,
  biExtremeHillsEdge,
  biJungle,
  biJungleHills,
  biJungleEdge,
  biDeepOcean,
  biStoneBeach,
  biColdBeach,
  biBirchForest,
  biBirchForestHills,
  biRoofedForest,
  biColdTaiga,
  biColdTaigaHills,
  biMegaTaiga,
  biMegaTaigaHills,
  biExtremeHillsPlus,
  biSavanna,
  biSavannaPlateau,
  biMesa,
  biMesaPlateauF,
  biMesaPlateau,
  biSunflowerPlains,
  biDesertM,
  biExtremeHillsM,
  biFlowerForest,
  biTaigaM,
  biSwamplandM,
  biIcePlainsSpikes,
  biJungleM,
  biJungleEdgeM,
  biBirchForestM,
  biBirchForestHillsM,
  biRoofedForestM,
  biColdTaigaM,
  biMegaSpruceTaiga,
  biMegaSpruceTaigaHills,
  biExtremeHillsPlusM,
  biSavannaM,
  biSavannaPlateauM,
  biMesaBryce,
  biMesaPlateauFM,
  biMesaPlateauM,
  biLava,
  biFlowingRiver,

  // ! Num Biomes should be the last element in the list since it represents the total number of biomes
  NUM_BIOMES,
};

namespace BiomeHelper
{
  inline bool isLiquidBiome(uint8_t b)
  {
    return b == (uint8_t)BIOME::biOcean ||
           b == (uint8_t)BIOME::biRiver ||
           b == (uint8_t)BIOME::biFlowingRiver ||
           b == (uint8_t)BIOME::biSwampland ||
           b == (uint8_t)BIOME::biFrozenRiver ||
           b == (uint8_t)BIOME::biFrozenOcean ||
           b == (uint8_t)BIOME::biLava;
  }
};

class Biome {
  public:
    float baseHeight;
    float amps[3][2];
    unsigned int color;
    uint8_t index;
};

const BIOME BIOMES_TEMPERATURE_HUMIDITY[] {
		//       0                1                2                      3                      4                      5                      6                7                8                9                10                     11                    12                     13                     14                    15
		/*  0 */ BIOME::biTundra, BIOME::biTundra, BIOME::biTundra,       BIOME::biTundra,       BIOME::biPlains,       BIOME::biPlains,       BIOME::biPlains, BIOME::biPlains, BIOME::biDesert, BIOME::biDesert, BIOME::biDesert,       BIOME::biDesert,      BIOME::biDesert,       BIOME::biDesert,       BIOME::biDesert,      BIOME::biDesert,
		/*  1 */ BIOME::biTundra, BIOME::biTundra, BIOME::biTundra,       BIOME::biTundra,       BIOME::biPlains,       BIOME::biPlains,       BIOME::biPlains, BIOME::biPlains, BIOME::biPlains, BIOME::biDesert, BIOME::biDesert,       BIOME::biDesert,      BIOME::biDesert,       BIOME::biDesert,       BIOME::biDesert,      BIOME::biDesert,
		/*  2 */ BIOME::biTundra, BIOME::biTundra, BIOME::biTundra,       BIOME::biTundra,       BIOME::biPlains,       BIOME::biExtremeHills, BIOME::biPlains, BIOME::biPlains, BIOME::biPlains, BIOME::biDesert, BIOME::biDesert,       BIOME::biDesert,      BIOME::biDesertHills,  BIOME::biDesertHills,  BIOME::biDesert,      BIOME::biDesert,
		/*  3 */ BIOME::biTundra, BIOME::biTundra, BIOME::biTundra,       BIOME::biTundra,       BIOME::biExtremeHills, BIOME::biExtremeHills, BIOME::biPlains, BIOME::biPlains, BIOME::biPlains, BIOME::biPlains, BIOME::biDesert,       BIOME::biDesert,      BIOME::biDesertHills,  BIOME::biDesertHills,  BIOME::biDesert,      BIOME::biDesert,
		/*  4 */ BIOME::biTundra, BIOME::biTundra, BIOME::biIceMountains, BIOME::biIceMountains, BIOME::biExtremeHills, BIOME::biExtremeHills, BIOME::biPlains, BIOME::biPlains, BIOME::biPlains, BIOME::biPlains, BIOME::biForestHills,  BIOME::biDesert, BIOME::biDesert, BIOME::biDesert, BIOME::biDesertHills, BIOME::biDesert,
		/*  5 */ BIOME::biTundra, BIOME::biTundra, BIOME::biIceMountains, BIOME::biIceMountains, BIOME::biExtremeHills, BIOME::biExtremeHills, BIOME::biPlains, BIOME::biPlains, BIOME::biPlains, BIOME::biPlains, BIOME::biForestHills,  BIOME::biForestHills, BIOME::biDesert, BIOME::biDesert, BIOME::biDesertHills, BIOME::biDesert,
		/*  6 */ BIOME::biTundra, BIOME::biTundra, BIOME::biIceMountains, BIOME::biIceMountains, BIOME::biForestHills,  BIOME::biForestHills,  BIOME::biForest, BIOME::biForest, BIOME::biForest, BIOME::biForest, BIOME::biForest,       BIOME::biForestHills, BIOME::biExtremeHills, BIOME::biDesert, BIOME::biDesertHills,      BIOME::biDesert,
		/*  7 */ BIOME::biTundra, BIOME::biTundra, BIOME::biIceMountains, BIOME::biIceMountains, BIOME::biForestHills,  BIOME::biForestHills,  BIOME::biForest, BIOME::biForest, BIOME::biForest, BIOME::biForest, BIOME::biForest,       BIOME::biForestHills, BIOME::biExtremeHills, BIOME::biExtremeHills, BIOME::biDesertHills,      BIOME::biDesert,
		/*  8 */ BIOME::biTundra, BIOME::biTundra, BIOME::biTaiga,        BIOME::biTaiga,        BIOME::biForest,       BIOME::biForest,       BIOME::biForest, BIOME::biForest, BIOME::biForest, BIOME::biForest, BIOME::biForest,       BIOME::biForestHills, BIOME::biExtremeHills, BIOME::biExtremeHills, BIOME::biPlains,      BIOME::biDesert,
		/*  9 */ BIOME::biTundra, BIOME::biTundra, BIOME::biTaiga,        BIOME::biTaiga,        BIOME::biForest,       BIOME::biForest,       BIOME::biForest, BIOME::biForest, BIOME::biForest, BIOME::biForest, BIOME::biForest,       BIOME::biForestHills, BIOME::biExtremeHills, BIOME::biExtremeHills, BIOME::biPlains,      BIOME::biPlains,
		/* 10 */ BIOME::biTaiga,  BIOME::biTaiga,  BIOME::biTaiga,        BIOME::biIceMountains, BIOME::biForestHills,  BIOME::biForestHills,  BIOME::biForest, BIOME::biForest, BIOME::biForest, BIOME::biForest, BIOME::biJungle,       BIOME::biJungle,      BIOME::biSwampland,    BIOME::biSwampland,    BIOME::biSwampland,   BIOME::biSwampland,
		/* 11 */ BIOME::biTaiga,  BIOME::biTaiga,  BIOME::biIceMountains, BIOME::biIceMountains, BIOME::biExtremeHills, BIOME::biForestHills,  BIOME::biForest, BIOME::biForest, BIOME::biForest, BIOME::biForest, BIOME::biJungle,       BIOME::biJungle,      BIOME::biSwampland,    BIOME::biSwampland,    BIOME::biSwampland,   BIOME::biSwampland,
		/* 12 */ BIOME::biTaiga,  BIOME::biTaiga,  BIOME::biIceMountains, BIOME::biIceMountains, BIOME::biExtremeHills, BIOME::biJungleHills,  BIOME::biJungle, BIOME::biJungle, BIOME::biJungle, BIOME::biJungle, BIOME::biJungle,       BIOME::biJungle,      BIOME::biSwampland,    BIOME::biSwampland,    BIOME::biSwampland,   BIOME::biSwampland,
		/* 13 */ BIOME::biTaiga,  BIOME::biTaiga,  BIOME::biTaiga,        BIOME::biIceMountains, BIOME::biJungleHills,  BIOME::biJungleHills,  BIOME::biJungle, BIOME::biJungle, BIOME::biJungle, BIOME::biJungle, BIOME::biJungle,       BIOME::biJungle,      BIOME::biSwampland,    BIOME::biSwampland,    BIOME::biSwampland,   BIOME::biSwampland,
		/* 14 */ BIOME::biTaiga,  BIOME::biTaiga,  BIOME::biTaiga,        BIOME::biTaiga,        BIOME::biJungle,       BIOME::biJungle,       BIOME::biJungle, BIOME::biJungle, BIOME::biJungle, BIOME::biJungle, BIOME::biJungle,       BIOME::biJungle,      BIOME::biSwampland,    BIOME::biSwampland,    BIOME::biSwampland,   BIOME::biSwampland,
		/* 15 */ BIOME::biTaiga,  BIOME::biTaiga,  BIOME::biTaiga,        BIOME::biTaiga,        BIOME::biJungle,       BIOME::biJungle,       BIOME::biJungle, BIOME::biJungle, BIOME::biJungle, BIOME::biJungle, BIOME::biJungle,       BIOME::biJungle,      BIOME::biSwampland,    BIOME::biSwampland,    BIOME::biSwampland,   BIOME::biSwampland,
	};


enum class MATERIAL : uint8_t {
  GRASS,
  DIRT,
  ROCK,
  STONE,

  // ! Num Materials should be the last element in the enum since it is used to determine the number of materials
  NUM_MATERIALS,
};

#endif // _BIOMES_H_