// example graph data for 6 procs
struct exGraphData {
static int exVertglbtabRaw[7];

static int exVertloctabRaw[69];
static int exVertwgtcmpRaw[69];
static int exVertwgtvisRaw[69];

static int exEdgeloctabRaw[220];
static int exEdgewgtcmpRaw[220];
static int exEdgewgtvisRaw[220];
};

int exGraphData::exVertglbtabRaw[7] = {   0,   10,   21,   32,   42,   53,   63};

int exGraphData::exVertloctabRaw[69] = {   0,    2,    5,    8,   11,   14,   17,   20,   23,   25,   28,    0,    4,    8,   12,   16,   20,   24,   28,   31,   34,   38,   42,    0,    4,    8,   12,   16,   20,   23,   26,   30,   34,   38,   42,    0,    4,    8,   12,   15,   18,   22,   26,   30,   34,   38,    0,    4,    8,   11,   14,   18,   22,   26,   30,   34,   38,   42,    0,    3,    5,    8,   11,   14,   17,   20,   23,   26,   28};
int exGraphData::exVertwgtcmpRaw[69] = {1100, 1101, 1102, 1103, 1104, 1105, 1106, 1107, 1108, 1109, 1110, 1110, 1111, 1112, 1113, 1114, 1115, 1116, 1117, 1118, 1119, 1120, 1121, 1120, 1121, 1122, 1123, 1124, 1125, 1126, 1127, 1128, 1129, 1130, 1131, 1130, 1131, 1132, 1133, 1134, 1135, 1136, 1137, 1138, 1139, 1140, 1140, 1141, 1142, 1143, 1144, 1145, 1146, 1147, 1148, 1149, 1150, 1151, 1150, 1151, 1152, 1153, 1154, 1155, 1156, 1157, 1158, 1159, 1160};
int exGraphData::exVertwgtvisRaw[69] = {1200, 1201, 1202, 1203, 1204, 1205, 1206, 1207, 1208, 1209, 1210, 1210, 1211, 1212, 1213, 1214, 1215, 1216, 1217, 1218, 1219, 1220, 1221, 1220, 1221, 1222, 1223, 1224, 1225, 1226, 1227, 1228, 1229, 1230, 1231, 1230, 1231, 1232, 1233, 1234, 1235, 1236, 1237, 1238, 1239, 1240, 1240, 1241, 1242, 1243, 1244, 1245, 1246, 1247, 1248, 1249, 1250, 1251, 1250, 1251, 1252, 1253, 1254, 1255, 1256, 1257, 1258, 1259, 1260};

int exGraphData::exEdgeloctabRaw[220] = {   1,    9,    0,    2,   10,    1,    3,   11,    2,    4,   12,    3,    5,   13,    4,    6,   14,    5,    7,   15,    6,    8,   16,    7,   17,    0,   10,   18,    1,    9,   11,   19,    2,   10,   12,   20,    3,   11,   13,   21,    4,   12,   14,   22,    5,   13,   15,   23,    6,   14,   16,   24,    7,   15,   17,   25,    8,   16,   26,    9,   19,   27,   10,   18,   20,   28,   11,   19,   21,   29,   12,   20,   22,   30,   13,   21,   23,   31,   14,   22,   24,   32,   15,   23,   25,   33,   16,   24,   26,   34,   17,   25,   35,   18,   28,   36,   19,   27,   29,   37,   20,   28,   30,   38,   21,   29,   31,   39,   22,   30,   32,   40,   23,   31,   33,   41,   24,   32,   34,   42,   25,   33,   35,   43,   26,   34,   44,   27,   37,   45,   28,   36,   38,   46,   29,   37,   39,   47,   30,   38,   40,   48,   31,   39,   41,   49,   32,   40,   42,   50,   33,   41,   43,   51,   34,   42,   44,   52,   35,   43,   53,   36,   46,   54,   37,   45,   47,   55,   38,   46,   48,   56,   39,   47,   49,   57,   40,   48,   50,   58,   41,   49,   51,   59,   42,   50,   52,   60,   43,   51,   53,   61,   44,   52,   62,   45,   55,   46,   54,   56,   47,   55,   57,   48,   56,   58,   49,   57,   59,   50,   58,   60,   51,   59,   61,   52,   60,   62,   53,   61};
int exGraphData::exEdgewgtcmpRaw[220] = {1100, 1101, 1100, 1102, 1103, 1102, 1104, 1105, 1104, 1106, 1107, 1106, 1108, 1109, 1108, 1110, 1111, 1110, 1112, 1113, 1112, 1114, 1115, 1114, 1116, 1101, 1117, 1118, 1103, 1117, 1119, 1120, 1105, 1119, 1121, 1122, 1107, 1121, 1123, 1124, 1109, 1123, 1125, 1126, 1111, 1125, 1127, 1128, 1113, 1127, 1129, 1130, 1115, 1129, 1131, 1132, 1116, 1131, 1133, 1118, 1134, 1135, 1120, 1134, 1136, 1137, 1122, 1136, 1138, 1139, 1124, 1138, 1140, 1141, 1126, 1140, 1142, 1143, 1128, 1142, 1144, 1145, 1130, 1144, 1146, 1147, 1132, 1146, 1148, 1149, 1133, 1148, 1150, 1135, 1151, 1152, 1137, 1151, 1153, 1154, 1139, 1153, 1155, 1156, 1141, 1155, 1157, 1158, 1143, 1157, 1159, 1160, 1145, 1159, 1161, 1162, 1147, 1161, 1163, 1164, 1149, 1163, 1165, 1166, 1150, 1165, 1167, 1152, 1168, 1169, 1154, 1168, 1170, 1171, 1156, 1170, 1172, 1173, 1158, 1172, 1174, 1175, 1160, 1174, 1176, 1177, 1162, 1176, 1178, 1179, 1164, 1178, 1180, 1181, 1166, 1180, 1182, 1183, 1167, 1182, 1184, 1169, 1185, 1186, 1171, 1185, 1187, 1188, 1173, 1187, 1189, 1190, 1175, 1189, 1191, 1192, 1177, 1191, 1193, 1194, 1179, 1193, 1195, 1196, 1181, 1195, 1197, 1198, 1183, 1197, 1199, 1200, 1184, 1199, 1201, 1186, 1202, 1188, 1202, 1203, 1190, 1203, 1204, 1192, 1204, 1205, 1194, 1205, 1206, 1196, 1206, 1207, 1198, 1207, 1208, 1200, 1208, 1209, 1201, 1209};
int exGraphData::exEdgewgtvisRaw[220] = {1200, 1201, 1200, 1202, 1203, 1202, 1204, 1205, 1204, 1206, 1207, 1206, 1208, 1209, 1208, 1210, 1211, 1210, 1212, 1213, 1212, 1214, 1215, 1214, 1216, 1201, 1217, 1218, 1203, 1217, 1219, 1220, 1205, 1219, 1221, 1222, 1207, 1221, 1223, 1224, 1209, 1223, 1225, 1226, 1211, 1225, 1227, 1228, 1213, 1227, 1229, 1230, 1215, 1229, 1231, 1232, 1216, 1231, 1233, 1218, 1234, 1235, 1220, 1234, 1236, 1237, 1222, 1236, 1238, 1239, 1224, 1238, 1240, 1241, 1226, 1240, 1242, 1243, 1228, 1242, 1244, 1245, 1230, 1244, 1246, 1247, 1232, 1246, 1248, 1249, 1233, 1248, 1250, 1235, 1251, 1252, 1237, 1251, 1253, 1254, 1239, 1253, 1255, 1256, 1241, 1255, 1257, 1258, 1243, 1257, 1259, 1260, 1245, 1259, 1261, 1262, 1247, 1261, 1263, 1264, 1249, 1263, 1265, 1266, 1250, 1265, 1267, 1252, 1268, 1269, 1254, 1268, 1270, 1271, 1256, 1270, 1272, 1273, 1258, 1272, 1274, 1275, 1260, 1274, 1276, 1277, 1262, 1276, 1278, 1279, 1264, 1278, 1280, 1281, 1266, 1280, 1282, 1283, 1267, 1282, 1284, 1269, 1285, 1286, 1271, 1285, 1287, 1288, 1273, 1287, 1289, 1290, 1275, 1289, 1291, 1292, 1277, 1291, 1293, 1294, 1279, 1293, 1295, 1296, 1281, 1295, 1297, 1298, 1283, 1297, 1299, 1300, 1284, 1299, 1301, 1286, 1302, 1288, 1302, 1303, 1290, 1303, 1304, 1292, 1304, 1305, 1294, 1305, 1306, 1296, 1306, 1307, 1298, 1307, 1308, 1300, 1308, 1309, 1301, 1309};
