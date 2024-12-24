import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Objects;

public class RandWalk {

    private static Matrix readMatrix(String filePath) throws IOException {
        ArrayList<ArrayList<Double>> result = new ArrayList<>();
        String[] lines = Files.readAllLines(Paths.get(filePath)).toArray(new String[0]);
        for (int i = 0; i < lines.length; i++) {
            result.add(new ArrayList<>());
            String[] values = lines[i].split("\t");
            for (String s : values) {
                double value = Double.parseDouble(s);
                result.get(i).add(value);
            }
        }
        return new Matrix(result);
    }

    // 归一化概率转移矩阵 pm
    public static Matrix normalizeMatrix(Matrix matrix) throws Exception {
        for (int i = 0; i < matrix.getRowCount(); i++) {
            double rowSum = 0.0;
            // 计算当前行的和
            for (int j = 0; j < matrix.getColCount(); j++) {
                rowSum += matrix.getMatrixValue(i, j);
            }
            // 如果该行和不为零则进行归一化
            if (rowSum != 0) {
                for (int j = 0; j < matrix.getColCount(); j++) {
                    matrix.setMatrixValue(i, j, matrix.getMatrixValue(i, j) / rowSum);
                }
            }
        }
        return matrix;
    }

    // 在随机游走图上游走，计算基因之间的相关性得分矩阵
    public static void walk(String gene_type, double arg) throws Exception {
        // 根据参数修改文件路径
        final ArrayList<String> filePaths = getPaths(gene_type);
        Matrix gm = readMatrix(filePaths.get(0)); // 原始矩阵
        Matrix pm = normalizeMatrix(readMatrix(filePaths.get(1))); // 概率转移矩阵，执行归一化
        Matrix pt = pm.transpose();
        // 计算基因相似度得分矩阵
        for (int k = 0; k < gm.getRowCount(); k++) {
            // RWR计算公式
            for (int i = 0; i < 15; i++) {
                Matrix newRow = Matrix.times(gm.getRow(k), pt);
                Matrix.dot(newRow, arg);
                newRow = Matrix.plus(newRow, Matrix.dot(Matrix.getEVector(pt.getColCount(), k), 1.0D - arg));
                Matrix.setRow(k, newRow, gm);
            }
        }
        // 基因相关性得分矩阵写入result0.x.mat文件
        final Path path = Paths.get("result/" + gene_type + "/result" + arg + ".mat");
        Calculator.createFileIfNotExists(path);
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < gm.getRowCount(); i++) {
            for (int j = 0; j < gm.getColCount(); j++)
                sb.append(gm.getMatrixValue(i, j)).append("\t");
            sb.append("\n");
            Files.write(path, sb.toString().getBytes(), StandardOpenOption.APPEND);
            sb = new StringBuilder();
        }
    }

    // 获取文件路径
    private static ArrayList<String> getPaths(String gene_type) {
        ArrayList<String> filePaths = new ArrayList<>();
        if (Objects.equals(gene_type, "Yeast")) {
            filePaths.add("temp/Yeast/netMatrix.txt");
            filePaths.add("temp/Yeast/pMatrix.txt");
        } else if (Objects.equals(gene_type, "Human")) {
            filePaths.add("temp/Human/netMatrix.txt");
            filePaths.add("temp/Human/pMatrix.txt");
        } else {
            filePaths.add("temp/Arabidopsis_thaliana/netMatrix.txt");
            filePaths.add("temp/Arabidopsis_thaliana/pMatrix.txt");
        }
        return filePaths;
    }
}
