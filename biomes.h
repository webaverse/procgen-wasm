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

// ? INFO : Biome types determine the terrain height function 
// TODO : Add more biome types 🙂
enum class BIOME : uint8_t {
  NULL_BIOME,

  DESERT_MOUNTAINS,
  FOREST_MOUNTAINS,
  ICE_MOUNTAINS,

  // Num Biomes should be the last element in the list since it represents the total number of biomes
  NUM_BIOMES,
};

// ? INFO : Liquid types determine the liquid height function
enum class LIQUID : uint8_t {
  NULL_LIQUID,

  OCEAN,
  RIVER,
  LAVA,
  FLOWING_RIVER,

  // Swamp

  // Num Liquids should be the last element in the list since it represents the total number of biomes
  NUM_LIQUIDS,
};

enum class MATERIAL : uint8_t {
  GRASS,
  DIRT,
  ROCK,
  STONE,

  // Num Materials should be the last element in the enum since it is used to determine the number of materials
  NUM_MATERIALS,
};

enum class INSTANCE : uint8_t {
  NULL_INSTANCE,

  TREE,
  BUSH,
  ROCK,
  STONE,
  GRASS,

  // Num Instances should be the last element in the list since it represents the total number of biomes
  NUM_INSTANCES,
};

enum class TREE : uint8_t {
  NULL_TREE,

  SHORT_TREE,
  MEDIUM_TREE,
  TALL_TREE,

  // Num Trees should be the last element in the list since it represents the total number of biomes
  NUM_TREES,
};

// class Biome {
//   public:
//     float baseHeight;
//     float amps[3][2];
//     unsigned int color;
//     uint8_t index;
// };
// const enum class BIOME : uint8_t {
//   biPlains,
//   biDesert,
//   biExtremeHills,
//   biForest,
//   biTaiga,
//   biNether,
//   biEnd,
//   biTundra,
//   biIceMountains,
//   biMushroomIsland,
//   biMushroomShore,
//   biBeach,
//   biDesertHills,
//   biForestHills,
//   biTaigaHills,
//   biExtremeHillsEdge,
//   biJungle,
//   biJungleHills,
//   biJungleEdge,
//   biStoneBeach,
//   biColdBeach,
//   biBirchForest,
//   biBirchForestHills,
//   biRoofedForest,
//   biColdTaiga,
//   biColdTaigaHills,
//   biMegaTaiga,
//   biMegaTaigaHills,
//   biExtremeHillsPlus,
//   biSavanna,
//   biSavannaPlateau,
//   biMesa,
//   biMesaPlateauF,
//   biMesaPlateau,
//   biSunflowerPlains,
//   biDesertM,
//   biExtremeHillsM,
//   biFlowerForest,
//   biTaigaM,
//   biSwamplandM,
//   biIcePlainsSpikes,
//   biJungleM,
//   biJungleEdgeM,
//   biBirchForestM,
//   biBirchForestHillsM,
//   biRoofedForestM,
//   biColdTaigaM,
//   biMegaSpruceTaiga,
//   biMegaSpruceTaigaHills,
//   biExtremeHillsPlusM,
//   biSavannaM,
//   biSavannaPlateauM,
//   biMesaBryce,
//   biMesaPlateauFM,
//   biMesaPlateauM,

//   NULL,

//   // Num Biomes should be the last element in the list since it represents the total number of biomes
//   NUM_BIOMES,
// };

// const BIOME BIOMES_TEMPERATURE_HUMIDITY[] {
// 		//       0                1                2                      3                      4                      5                      6                7                8                9                10                     11                    12                     13                     14                    15
// 		/*  0 */ BIOME::biTundra, BIOME::biTundra, BIOME::biTundra,       BIOME::biTundra,       BIOME::biPlains,       BIOME::biPlains,       BIOME::biPlains, BIOME::biPlains, BIOME::biDesert, BIOME::biDesert, BIOME::biDesert,       BIOME::biDesert,      BIOME::biDesert,       BIOME::biDesert,       BIOME::biDesert,      BIOME::biDesert,
// 		/*  1 */ BIOME::biTundra, BIOME::biTundra, BIOME::biTundra,       BIOME::biTundra,       BIOME::biPlains,       BIOME::biPlains,       BIOME::biPlains, BIOME::biPlains, BIOME::biPlains, BIOME::biDesert, BIOME::biDesert,       BIOME::biDesert,      BIOME::biDesert,       BIOME::biDesert,       BIOME::biDesert,      BIOME::biDesert,
// 		/*  2 */ BIOME::biTundra, BIOME::biTundra, BIOME::biTundra,       BIOME::biTundra,       BIOME::biPlains,       BIOME::biExtremeHills, BIOME::biPlains, BIOME::biPlains, BIOME::biPlains, BIOME::biDesert, BIOME::biDesert,       BIOME::biDesert,      BIOME::biDesertHills,  BIOME::biDesertHills,  BIOME::biDesert,      BIOME::biDesert,
// 		/*  3 */ BIOME::biTundra, BIOME::biTundra, BIOME::biTundra,       BIOME::biTundra,       BIOME::biExtremeHills, BIOME::biExtremeHills, BIOME::biPlains, BIOME::biPlains, BIOME::biPlains, BIOME::biPlains, BIOME::biDesert,       BIOME::biDesert,      BIOME::biDesertHills,  BIOME::biDesertHills,  BIOME::biDesert,      BIOME::biDesert,
// 		/*  4 */ BIOME::biTundra, BIOME::biTundra, BIOME::biIceMountains, BIOME::biIceMountains, BIOME::biExtremeHills, BIOME::biExtremeHills, BIOME::biPlains, BIOME::biPlains, BIOME::biPlains, BIOME::biPlains, BIOME::biForestHills,  BIOME::biDesert, BIOME::biDesert, BIOME::biDesert, BIOME::biDesertHills, BIOME::biDesert,
// 		/*  5 */ BIOME::biTundra, BIOME::biTundra, BIOME::biIceMountains, BIOME::biIceMountains, BIOME::biExtremeHills, BIOME::biExtremeHills, BIOME::biPlains, BIOME::biPlains, BIOME::biPlains, BIOME::biPlains, BIOME::biForestHills,  BIOME::biForestHills, BIOME::biDesert, BIOME::biDesert, BIOME::biDesertHills, BIOME::biDesert,
// 		/*  6 */ BIOME::biTundra, BIOME::biTundra, BIOME::biIceMountains, BIOME::biIceMountains, BIOME::biForestHills,  BIOME::biForestHills,  BIOME::biForest, BIOME::biForest, BIOME::biForest, BIOME::biForest, BIOME::biForest,       BIOME::biForestHills, BIOME::biExtremeHills, BIOME::biDesert, BIOME::biDesertHills,      BIOME::biDesert,
// 		/*  7 */ BIOME::biTundra, BIOME::biTundra, BIOME::biIceMountains, BIOME::biIceMountains, BIOME::biForestHills,  BIOME::biForestHills,  BIOME::biForest, BIOME::biForest, BIOME::biForest, BIOME::biForest, BIOME::biForest,       BIOME::biForestHills, BIOME::biExtremeHills, BIOME::biExtremeHills, BIOME::biDesertHills,      BIOME::biDesert,
// 		/*  8 */ BIOME::biTundra, BIOME::biTundra, BIOME::biTaiga,        BIOME::biTaiga,        BIOME::biForest,       BIOME::biForest,       BIOME::biForest, BIOME::biForest, BIOME::biForest, BIOME::biForest, BIOME::biForest,       BIOME::biForestHills, BIOME::biExtremeHills, BIOME::biExtremeHills, BIOME::biPlains,      BIOME::biDesert,
// 		/*  9 */ BIOME::biTundra, BIOME::biTundra, BIOME::biTaiga,        BIOME::biTaiga,        BIOME::biForest,       BIOME::biForest,       BIOME::biForest, BIOME::biForest, BIOME::biForest, BIOME::biForest, BIOME::biForest,       BIOME::biForestHills, BIOME::biExtremeHills, BIOME::biExtremeHills, BIOME::biPlains,      BIOME::biPlains,
// 		/* 10 */ BIOME::biTaiga,  BIOME::biTaiga,  BIOME::biTaiga,        BIOME::biIceMountains, BIOME::biForestHills,  BIOME::biForestHills,  BIOME::biForest, BIOME::biForest, BIOME::biForest, BIOME::biForest, BIOME::biJungle,       BIOME::biJungle,      BIOME::biSwampland,    BIOME::biSwampland,    BIOME::biSwampland,   BIOME::biSwampland,
// 		/* 11 */ BIOME::biTaiga,  BIOME::biTaiga,  BIOME::biIceMountains, BIOME::biIceMountains, BIOME::biExtremeHills, BIOME::biForestHills,  BIOME::biForest, BIOME::biForest, BIOME::biForest, BIOME::biForest, BIOME::biJungle,       BIOME::biJungle,      BIOME::biSwampland,    BIOME::biSwampland,    BIOME::biSwampland,   BIOME::biSwampland,
// 		/* 12 */ BIOME::biTaiga,  BIOME::biTaiga,  BIOME::biIceMountains, BIOME::biIceMountains, BIOME::biExtremeHills, BIOME::biJungleHills,  BIOME::biJungle, BIOME::biJungle, BIOME::biJungle, BIOME::biJungle, BIOME::biJungle,       BIOME::biJungle,      BIOME::biSwampland,    BIOME::biSwampland,    BIOME::biSwampland,   BIOME::biSwampland,
// 		/* 13 */ BIOME::biTaiga,  BIOME::biTaiga,  BIOME::biTaiga,        BIOME::biIceMountains, BIOME::biJungleHills,  BIOME::biJungleHills,  BIOME::biJungle, BIOME::biJungle, BIOME::biJungle, BIOME::biJungle, BIOME::biJungle,       BIOME::biJungle,      BIOME::biSwampland,    BIOME::biSwampland,    BIOME::biSwampland,   BIOME::biSwampland,
// 		/* 14 */ BIOME::biTaiga,  BIOME::biTaiga,  BIOME::biTaiga,        BIOME::biTaiga,        BIOME::biJungle,       BIOME::biJungle,       BIOME::biJungle, BIOME::biJungle, BIOME::biJungle, BIOME::biJungle, BIOME::biJungle,       BIOME::biJungle,      BIOME::biSwampland,    BIOME::biSwampland,    BIOME::biSwampland,   BIOME::biSwampland,
// 		/* 15 */ BIOME::biTaiga,  BIOME::biTaiga,  BIOME::biTaiga,        BIOME::biTaiga,        BIOME::biJungle,       BIOME::biJungle,       BIOME::biJungle, BIOME::biJungle, BIOME::biJungle, BIOME::biJungle, BIOME::biJungle,       BIOME::biJungle,      BIOME::biSwampland,    BIOME::biSwampland,    BIOME::biSwampland,   BIOME::biSwampland,
// 	};

#endif // _BIOMES_H_