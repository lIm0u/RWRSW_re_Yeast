import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.OpenOption;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.nio.file.attribute.FileAttribute;
import java.util.ArrayList;

public class RandWalk {
    // 读取概率转移矩阵P
    private static Matrix readPMatrix() throws IOException {
        ArrayList<ArrayList<Double>> result = new ArrayList<>();
        String[] lines = Files.readAllLines(Paths.get("temp/pMatrix.txt")).toArray(new String[0]);
        for (int i = 0; i < lines.length; i++) {
            result.add(new ArrayList<>(4172));
            String[] values = lines[i].split("\t");
            for (int j = 0; j < values.length; j++) {
                double value = Double.parseDouble(values[j]);
                result.get(i).add(value);
            }
        }
//        // 归一化处理
//        for(int i=0;i<lines.length;i++){
//            double sum = 0.0D;
//            int j;
//            for(j=0;j<lines.length;j++)
//                sum += ((result.get(j)).get(i)).doubleValue() * ((result.get(j)).get(i)).doubleValue();
//            for(j=0;j<lines.length;j++){
//                if(((result.get(j)).get(i)).doubleValue() != 0.0D)
//                    result.get(j).set(i,Double.valueOf(((result.get(j)).get(i)).doubleValue() / sum));
//            }
//        }
        return new Matrix(result);
    }

    // 读取原始网络图矩阵
    private static Matrix readGMatrix() throws IOException {
        ArrayList<ArrayList<Double>> result = new ArrayList<>();
        String[] lines = Files.readAllLines(Paths.get("temp/netMatrix.txt", new String[0])).toArray(new String[0]);
        for (int i = 0; i < lines.length; i++) {
            result.add(new ArrayList<>(4172));
            String[] values = lines[i].split("\t");
            for (int j = 0; j < values.length; j++) {
                double value = Double.parseDouble(values[j]);
                result.get(i).add(value);
            }
        }
        return new Matrix(result);
    }

    // 在随机游走图上游走，计算基因之间的相关性得分矩阵
    public static void walk(double arg) throws Exception {
        Matrix gm = readGMatrix(); // 原始矩阵
        Matrix pm = readPMatrix(); // 归一化的概率转移矩阵
        double c = arg; // 阿尔法参数

        long s = System.currentTimeMillis();
        // 计算基因相似度得分矩阵
        for (int k = 0; k < gm.getRowCount(); k++) {
            long start = System.currentTimeMillis();
            // RWR计算公式
            for (int i = 0; i < 15; i++) {
                Matrix newRow = Matrix.times(gm.getRow(k), pm);
                newRow = Matrix.dot(newRow, c);
                newRow = Matrix.plus(newRow, Matrix.dot(Matrix.getEVector(pm.getColCount(), k), 1.0D - c));
                Matrix.setRow(k, newRow, gm);
            }
            long end = System.currentTimeMillis();
            System.out.println("第" + k + "次，耗时：" + (end - start) + "总耗时：" + (end - s));
        }

        // 基因相关性得分矩阵写入result0.x.mat文件
        if (!Files.exists(Paths.get("result/result" + arg + ".mat"), new java.nio.file.LinkOption[0]))
            Files.createFile(Paths.get("result/result" + arg + ".mat"), new FileAttribute[0]);
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < 4172; i++) {
            for (int j = 0; j < 4172; j++)
                sb.append(gm.getMatrixValue(i, j) + "\t");
            sb.append("\n");
            Files.write(Paths.get("result/result" + arg + ".mat"), sb.toString().getBytes(), new OpenOption[]{StandardOpenOption.APPEND});
            sb = new StringBuilder();
        }
    }
}
