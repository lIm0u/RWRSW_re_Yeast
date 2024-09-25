import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.OpenOption;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.List;

public class Main {
    public static void main(String[] args) throws Exception {
        // 1、随机游走的相关性得分矩阵
//        Calculator.getRWRMatrix(0.9D);
//        System.out.println("矩阵已算完！");
//        System.out.println("======================================");

        // 2、计算基因相似度
        String input = "data/lfcpair.txt";
        final String output = "result/similarityResult.txt";
        Calculator.data_initializer();
        System.out.println("基因相似度计算中...");
        final String[] lines = Files.readAllLines(Paths.get(input)).toArray(new String[0]);
        int thCnt = 4;
        List<Thread> threads = new ArrayList<>();
        for (int i = 0; i < thCnt; i++) {
            final int startIndex = i;
            final int finalThCnt = thCnt;
            Thread th = new Thread(() -> {
                long sss = System.currentTimeMillis();
                for (int j = startIndex; j < lines.length; j += finalThCnt) {
                    String g1 = lines[j].split(":")[0];
                    String g2 = lines[j].split(":")[1];
                    double sim = Calculator.get_similarity_of_genes(g1, g2);
                    System.out.print("Thread:" + (startIndex + 1) + "\t" + g1 + "\t" + g2 + "\t" + sim);
                    System.out.println("\tTime Use:" + (System.currentTimeMillis() - sss));
                    if (!Files.exists(Paths.get(output))) {
                        try {
                            Files.createFile(Paths.get(output));
                        } catch (IOException e) {
                            throw new RuntimeException(e);
                        }
                    }
                    try {
                        Files.write(Paths.get(output), (g1 + "\t" + g2 + "\t" + sim + "\n").getBytes(), StandardOpenOption.APPEND);
                    } catch (IOException e) {
                        throw new RuntimeException(e);
                    }
                }
            });
            th.start();
            threads.add(th);
        }
        for (Thread thread : threads) {
            thread.join();
        }
//         3、根据ec分组计算各组的LFC得分
//        Calculator.Evaluator.LFC();
        System.out.println("LFC得分已算完！");
        System.out.println("Successfully!");
    }
}
