/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package processamento;

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import static java.lang.Math.sqrt;
import static java.lang.Math.pow;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.Stack;
import javax.imageio.ImageIO;

/**
 *
 * @author Álef e Ketson
 * @link https://github.com/alef-sena/PDI.Redness
 */
public class ProcessamentoImagem {
    
    // 
    public static BufferedImage rednessDetection1(BufferedImage img) {
        
        // Carregar a imagem
        BufferedImage res = new BufferedImage(img.getWidth(), img.getHeight(), BufferedImage.TYPE_BYTE_GRAY);

        double maiorCol = -Double.MAX_VALUE;
        double menorCol = Double.MAX_VALUE;
        
        for (int i = 0; i < img.getWidth(); i++) {
            for (int j = 0; j < img.getHeight(); j++) {
                Color pixelColor = new Color(img.getRGB(i, j));
                

                int canalVermelho = pixelColor.getRed();
                int canalVerde = pixelColor.getGreen();
                int canalAzul = pixelColor.getBlue();

                double intensidadeVermelha = canalVermelho - 0.5 * canalVerde - 0.5 * canalAzul;
                double intensidadeNormalizada = (intensidadeVermelha - (-255.0)) / (255.0 - (-255.0));

                if (intensidadeNormalizada > 1.0) {
                    intensidadeNormalizada = 1.0;
                } else if (intensidadeNormalizada < 0.0) {
                    intensidadeNormalizada = 0.0;
                }
                
                int corNorm = (int) (intensidadeNormalizada * 255);

                Color novoPixel = new Color(corNorm, corNorm, corNorm);
                res.setRGB(i, j, novoPixel.getRGB());

                if (intensidadeNormalizada > maiorCol) maiorCol = intensidadeNormalizada;
                if (intensidadeNormalizada < menorCol) menorCol = intensidadeNormalizada;
            }
        }

        for (int i = 0; i < img.getWidth(); i++) {
            for (int j = 0; j < img.getHeight(); j++) {
                int corNorm = res.getRGB(i, j);
                double intensidade = (corNorm >> 16) & 0xFF;
                double cor_norm = normalizacao(intensidade, menorCol * 255, maiorCol * 255, 0, 255);
                //System.out.println(cor_norm);
                Color atual = new Color((int) cor_norm, (int) cor_norm, (int) cor_norm);
                res.setRGB(i, j, atual.getRGB());
            }
        }

        return res;
    }
    
    // 
    public static BufferedImage rednessDetection2(BufferedImage originalImage) {
        
        int width = originalImage.getWidth();
        int height = originalImage.getHeight();

        BufferedImage res = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);

        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {

                int rgb = originalImage.getRGB(i, j);
                int red = (rgb >> 16) & 0xFF;

                int grayValue = red;

                int grayPixel = (grayValue << 16) | (grayValue << 8) | grayValue;

                res.setRGB(i, j, grayPixel);
            }
        }

        return res;
    }
    
    //
    public static BufferedImage rednessFusion(BufferedImage img) {
        
        BufferedImage res = new BufferedImage(img.getWidth(), img.getHeight(), BufferedImage.TYPE_BYTE_GRAY);
                
        for (int i = 0; i < img.getWidth(); i++) {
            for (int j = 0; j < img.getHeight(); j++) {
                
                // Algoritmo 1 - Transformar a imagem para níveis de cinza
                double intensidadeNormalizada;
                Color pixelColor = new Color(img.getRGB(i, j));

                int canalVermelho = pixelColor.getRed();
                int canalVerde = pixelColor.getGreen();
                int canalAzul = pixelColor.getBlue();

                double intensidadeVermelha = canalVermelho - 0.5 * canalVerde - 0.5 * canalAzul;
                intensidadeNormalizada = (intensidadeVermelha - (-255.0)) / (255.0 - (-255.0));

                if (intensidadeNormalizada > 1.0) {
                    intensidadeNormalizada = 1.0;
                } else if (intensidadeNormalizada < 0.0) {
                    intensidadeNormalizada = 0.0;
                }
                
                // Algoritmo 2 - Transformar a imagem para níveis de cinza
                int rgb = img.getRGB(i, j);
                int red = (rgb >> 16) & 0xFF;
                
                int grayValue = red;
                
                int grayPixel = (grayValue << 16) | (grayValue << 8) | grayValue;
                
                double pixelNormalized = (double) grayPixel / 255.0;
                
                res.setRGB(i, j, (int) ((intensidadeNormalizada > pixelNormalized ? intensidadeNormalizada : pixelNormalized) * 255));
            }
        }
        
        return res;
    }

    //
    public static BufferedImage finalDetection(BufferedImage img){
        
        BufferedImage initialRedness = new BufferedImage(img.getWidth(), img.getHeight(), BufferedImage.TYPE_BYTE_GRAY);
        BufferedImage finalImage = new BufferedImage(img.getWidth(), img.getHeight(), BufferedImage.TYPE_BYTE_GRAY);

        int width = img.getWidth();
        int height = img.getHeight();
        
        double maiorCol = -Double.MAX_VALUE;
        double menorCol = Double.MAX_VALUE;
        
        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                
                Color pixelColor = new Color(img.getRGB(i, j));

                int canalVermelho = pixelColor.getRed();
                int canalVerde = pixelColor.getGreen();
                int canalAzul = pixelColor.getBlue();

                double intensidadeVermelha = canalVermelho - 0.5 * canalVerde - 0.5 * canalAzul;
                double intensidadeNormalizada = (intensidadeVermelha - (-255.0)) / (255.0 - (-255.0));

                if (intensidadeNormalizada > 1.0) {
                    intensidadeNormalizada = 1.0;
                } else if (intensidadeNormalizada < 0.0) {
                    intensidadeNormalizada = 0.0;
                }

                int corNorm = (int) (intensidadeNormalizada * 255);

                Color novoPixel = new Color(corNorm, corNorm, corNorm);
                initialRedness.setRGB(i, j, novoPixel.getRGB());

                if (intensidadeNormalizada > maiorCol) maiorCol = intensidadeNormalizada;
                if (intensidadeNormalizada < menorCol) menorCol = intensidadeNormalizada;
            }
        }
        
        for (int a = 0; a < img.getWidth(); a++) {
            for (int b = 0; b < img.getHeight(); b++) {
                int normalizedColor = img.getRGB(a, b);
                double intensidade = (normalizedColor >> 16) & 0xFF;
                double cor_norm = normalizacao(intensidade, menorCol * 255, maiorCol * 255, 0, 255);
                //System.out.println(cor_norm);
                Color atual = new Color((int) cor_norm, (int) cor_norm, (int) cor_norm);
                //System.out.println("Normalizando: " + atual.getRGB());
                initialRedness.setRGB(a, b, atual.getRGB());
            }
        }
        
        // aqui eu tenho initialRedness nivel de cinza
        
                // Binarizacao
        //Armazenando o tamanho da imagem
        int Largura = initialRedness.getWidth();    // Largura da imagem
        int Altura = initialRedness.getHeight();    // Altura da imagem
        int col, lin, c, cinza;
        double totalPixel = (double) Largura*Altura;
        double [] proba = new double[256];
        int [] histogram = new int[256];

        // Inicializacao variaveis
        int k, uiLimiar;        
        // inicializacao do Histograma
        for(c = 0; c < 256; c++) histogram[c]= 0;

        // Passo 1: calculo do Histograma  
        for( lin = 0; lin < Altura; lin++) {
            for( col = 0; col < Largura; col++) {
                Color x = new Color(initialRedness.getRGB(col,lin)); // leitura dos canais RGB de cada pixel
                // caso a imagem for colorida, gera uma imagem em ní­veis de cinza = (R+G+B)/3
                cinza = (int)((x.getGreen() + x.getRed() + x.getBlue())/3);		
                histogram[cinza]++; // Atualiza o histograma dos nÃ­veis de cinza (entre 0 e 255)
            }
        }

        // Alocacao das Matrizes
        double fSigmaBMax, fMiTotal;
        double [] fOmega = new double[256], fMi = new double[256] , fSigmaB = new double[256];

        //Passo 2: calculo das probabilidades a priori
        for (c = 0; c < 256; c++) {
                proba[c] = (double) ((histogram[c])/(double)(totalPixel));
                fOmega[c] = fMi[c] = 0.0;
        }  

        for (k = 0; k < 256; k++){
            for (c = 0; c < k; c++)   fOmega[k] += proba[c];
        }

        for (k = 0; k < 256; k++){
            for (c = 0; c < k; c++)   fMi[k] += (c + 1) * proba[c];
        }                    

        fMiTotal = fSigmaBMax = 0.0;  
        uiLimiar = 128; //inicializacao do valor de limiar  de Otsu

        for (c = 0; c < 256; c++)       fMiTotal += (c + 1) * proba[c];

        // calculo de fSigmaBMax para 1o ni­vel de cinza = 0 
        if ((fOmega[0] * (1 - fOmega[0])) != 0.0) {
                // calculo do criterio de binarizacao da tecnica de Otsu 
                fSigmaBMax = (  (fMiTotal * fOmega[0] - fMi[0]) * (fMiTotal * fOmega[0] - fMi[0]) ) / (fOmega[0] * (1 - fOmega[0]));
                uiLimiar = 0;
        }

        // calculo de fSigmaBMax para os outros ni­veis de cinza (entre 1 e 255)
        for (k = 1; k < 256; k++) {
            if ((fOmega[k] * (1 - fOmega[k])) != 0.0){
                // Busca do valor de limiar = valor de ni­vel de cinza que corresponde ao maior valor de fSigmaBMax (criterio da tecnica de Otsu )
                fSigmaB[k] = (  (fMiTotal * fOmega[k] - fMi[k]) * (fMiTotal * fOmega[k] - fMi[k]) ) / (fOmega[k] * (1 - fOmega[k]));
                if (fSigmaB[k] > fSigmaBMax){
                    fSigmaBMax = fSigmaB[k];
                    uiLimiar = (int) k;
                }
            }
        }

        // Visualizacao do valor de limiar de Otsu  em modo Debug
        // System.out.println(uiLimiar);

        //Cria a  imagem  binarizada 
        //Aloca a Matriz
        int [][] pBufferbinario = new int[Altura][Largura]; //Cria um PONTEIRO para a  imagem  binarizada

        for( lin = 0; lin < Altura; lin++) {
            for( col = 0; col < Largura; col++) {

                Color x = new Color(initialRedness.getRGB(col,lin)); // leitura dos canais RGB de cada pixel
                
                int rgb = img.getRGB(col, lin);
                int red = (rgb >> 16) & 0xFF;
                int green = (rgb >> 8) & 0xFF;
                int blue = rgb & 0xFF;
                
                Color pixelColor = new Color(img.getRGB(col, lin));
                //Color pixelColor = new Color(img.getRGB(col, lin));

                // caso a imagem for colorida, gera uma imagem em ni­veis de cinza = (R+G+B)/3
                cinza = (int)((x.getGreen() + x.getRed() + x.getBlue())/3);

                if (cinza < uiLimiar){
                    pBufferbinario[lin][col] = 0;
                }
                else{
                    pixelColor = new Color(255, 255, 255);
                    pBufferbinario[lin][col] = 1;
                }  //multiplicado depois por 255
                //System.out.println(pixel.getRGB());
                finalImage.setRGB(col, lin, pixelColor.getRGB());
            }
        }
                
        return finalImage;
    }

    // 
    public static BufferedImage detectarPeleHumana(BufferedImage img){
        
        BufferedImage res = new BufferedImage(img.getWidth(), img.getHeight(), img.getType());
	
        //Para realizar a análise de forma eficiente, precisa-se de um vetor para armazenar o índice de vermelhidao para cada pixel
        double[][] vet_col = new double[img.getWidth()][img.getHeight()];
        double maior_col = -Double.MAX_VALUE, menor_col = Double.MAX_VALUE;     //Valores que armazenarão a maior e a menor intensidade da imagem
        
        //Busca pelo índice de vermelhidão segundo a fórmula correspondente à medição de cor dos tomates
	for (int i = 0; i < img.getWidth(); i++) {
            for (int j = 0; j < img.getHeight(); j++) {
                
                //Vetor para o armazenamento das informações provenientes da conversão de RGB para L*a*b*
                double[] Lab = new double[3];
                Lab = RGBtoLab(img, i, j);
                double L = Lab[0], a = Lab[1], b = Lab[2];
                
                //Obtenção do croma, característica do espaço de cores L*a*b*
                double chroma = sqrt(pow(a, 2) + pow(b, 2));
                
                //Obtendo o h
                double h = Math.atan(b/a);
                double H = Math.toDegrees(h);
                
                //Índice de vermelhidão dos tomates, segundo a referência utilizada
                double cirg2 = (180-H)/(L/chroma);
                
                //Encontrar a maior e a menor intensidades para a normalização
                if(cirg2 > maior_col) maior_col = cirg2;
                if(cirg2 < menor_col) menor_col = cirg2;
                vet_col[i][j] = cirg2;
            }
        }
        
        //Geração da nova imagem após normalização em níveis de cinza
        for (int i = 0; i < img.getWidth(); i++) {
            for (int j = 0; j < img.getHeight(); j++) {
                
                double cor_norm = normalizacao(vet_col[i][j], menor_col, maior_col, 0, 255);
                Color atual = new Color((int) cor_norm, (int) cor_norm, (int) cor_norm);
                res.setRGB(i, j, atual.getRGB());
                
            }
        }
        
        //A imagem após análise e medição de vermelhidão
        return res;
    } 
    
    // Normalizacao em níveis de cinza para a imagem
    private static double normalizacao(double hist, double min_hist, double max_hist, double min_index, double max_index){
        hist = Math.max(min_hist, Math.min(max_hist, hist));
        //A equação de normalização disposta no conteúdo teórico
        return (max_index - min_index)*(hist - min_hist)/(max_hist - min_hist) + min_index;
    }
    
    // 
    public static BufferedImage BinarizacaoOtsu(BufferedImage img) {
	
	 //Cria a  imagem  de sai­da  
	// Aloca a Matriz da imagem  de sai­da  
	BufferedImage res = new BufferedImage(img.getWidth(), img.getHeight(), img.getType());

        //Armazenando o tamanho da imagem
        int Largura = img.getWidth();    // Largura da imagem
        int Altura = img.getHeight();    // Altura da imagem
        int col, lin, i, cinza;
        double totalPixel = (double) Largura*Altura;
        double [] proba = new double[256];
        int [] histogram = new int[256];

        // Inicializacao variaveis
   	int k, uiLimiar;        
        // inicializacao do Histograma
        for(i = 0; i < 256; i++) histogram[i]= 0;
 
        // Passo 1: calculo do Histograma  
	for( lin = 0; lin < Altura; lin++) {
            for( col = 0; col < Largura; col++) {
                Color x = new Color(img.getRGB(col,lin)); // leitura dos canais RGB de cada pixel
		// caso a imagem for colorida, gera uma imagem em ní­veis de cinza = (R+G+B)/3
		cinza = (int)((x.getGreen() + x.getRed() + x.getBlue())/3);		
		histogram[cinza]++; // Atualiza o histograma dos nÃ­veis de cinza (entre 0 e 255)
            }
	}
        
        // Alocacao das Matrizes
        double fSigmaBMax, fMiTotal;
        double [] fOmega = new double[256], fMi = new double[256] , fSigmaB = new double[256];
  
        //Passo 2: calculo das probabilidades a priori
        for (i = 0; i < 256; i++) {
		proba[i] = (double) ((histogram[i])/(double)(totalPixel));
		fOmega[i] = fMi[i] = 0.0;
        }  
 
	for (k = 0; k < 256; k++){
            for (i = 0; i < k; i++)   fOmega[k] += proba[i];
        }

	for (k = 0; k < 256; k++){
            for (i = 0; i < k; i++)   fMi[k] += (i + 1) * proba[i];
        }                    
  
        fMiTotal = fSigmaBMax = 0.0;  
        uiLimiar = 128; //inicializacao do valor de limiar  de Otsu

        for (i = 0; i < 256; i++)       fMiTotal += (i + 1) * proba[i];
  
        // calculo de fSigmaBMax para 1o ni­vel de cinza = 0 
        if ((fOmega[0] * (1 - fOmega[0])) != 0.0) {
		// calculo do criterio de binarizacao da tecnica de Otsu 
		fSigmaBMax = (  (fMiTotal * fOmega[0] - fMi[0]) * (fMiTotal * fOmega[0] - fMi[0]) ) / (fOmega[0] * (1 - fOmega[0]));
		uiLimiar = 0;
        }
        
        // calculo de fSigmaBMax para os outros ni­veis de cinza (entre 1 e 255)
        for (k = 1; k < 256; k++) {
            if ((fOmega[k] * (1 - fOmega[k])) != 0.0){
                // Busca do valor de limiar = valor de ni­vel de cinza que corresponde ao maior valor de fSigmaBMax (criterio da tecnica de Otsu )
                fSigmaB[k] = (  (fMiTotal * fOmega[k] - fMi[k]) * (fMiTotal * fOmega[k] - fMi[k]) ) / (fOmega[k] * (1 - fOmega[k]));
		if (fSigmaB[k] > fSigmaBMax){
                    fSigmaBMax = fSigmaB[k];
                    uiLimiar = (int) k;
                }
            }
        }
    
        // Visualizacao do valor de limiar de Otsu  em modo Debug
        // System.out.println(uiLimiar);
    
        //Cria a  imagem  binarizada 
	//Aloca a Matriz
	int [][] pBufferbinario = new int[Altura][Largura]; //Cria um PONTEIRO para a  imagem  binarizada 

        for( lin = 0; lin < Altura; lin++) {
            for( col = 0; col < Largura; col++) {
                Color x = new Color(img.getRGB(col,lin)); // leitura dos canais RGB de cada pixel
                // caso a imagem for colorida, gera uma imagem em ni­veis de cinza = (R+G+B)/3
		cinza = (int)((x.getGreen() + x.getRed() + x.getBlue())/3);
			
		if (cinza < uiLimiar)   pBufferbinario[lin][col] = 0;
		else    pBufferbinario[lin][col] = 1;  //multiplicado depois por 255
            }
	}

        //Aqui Gera a  imagem binaria 
        for(lin = 0; lin < Altura; lin++) {
		for(col = 0; col < Largura; col++) {
                    int atual = pBufferbinario[lin][col]* 255; //multiplicacao da img binaria por 255
                    Color novo = new Color(atual, atual, atual);
                    res.setRGB(col,lin, novo.getRGB()); // gravacao do valor binarizado para cada pixel
		}
        }

    return res;

}
    
    // 
    public static BufferedImage EntropiaPun(BufferedImage img) {
	
	 //Cria a  imagem  de sai­da  
	// Aloca a Matriz da imagem  de sai­da  
	BufferedImage res = new BufferedImage(img.getWidth(), img.getHeight(), img.getType());

        //Armazenando o tamanho da imagem
        int Largura = img.getWidth();    // Largura da imagem
        int Altura = img.getHeight();    // Altura da imagem
        int col, lin, i, c, cinza;
        double totalPixel = (double) Largura*Altura;
        double [] proba = new double[256];
        double [] ni = new double[256];
        int [] histogram = new int[256];
        double HT, Ht, Pt;
        double MaxProb1, MaxProb2, Aux1, Aux2, Aux3, max;

        // Inicializacao variaveis
   	int uiLimiar = 0;        
        // inicializacao do Histograma
        for(i = 0; i < 256; i++) histogram[i]= 0;
 
        // Passo 1: calculo do Histograma  
	for( lin = 0; lin < Altura; lin++) {
            for( col = 0; col < Largura; col++) {
                Color x = new Color(img.getRGB(col,lin)); // leitura dos canais RGB de cada pixel
		// caso a imagem for colorida, gera uma imagem em ní­veis de cinza = (R+G+B)/3
		cinza = (int)((x.getGreen() + x.getRed() + x.getBlue())/3);		
		histogram[cinza]++; // Atualiza o histograma dos nÃ­veis de cinza (entre 0 e 255)
            }
	}
        
  
        //Passo 2: calculo das probabilidades a priori
        for (i = 0; i < 256; i++) {
		proba[i] = (double) ((histogram[i])/(double)(totalPixel));
		ni[i] = 0;
        }  
                    
        HT = 0;
        
        for(i = 0; i < 256; i++ ){
            if(proba[i] != 0){
                HT+= - proba[i] * Math.log(proba[i]);
            }
        }
        
        Pt = 0;
        Ht = 0;
        MaxProb1 = 0;
        
        for (i = 0; i < 256; i++){
            if(proba[i] != 0){
                
                //Entropia de Pixels Pretos
                Pt+= proba[i];
                Ht+= -proba[i] * Math.log(proba[i]);
                
                Aux3 = Ht/HT;
                
                // Max da Probabilidade de 1 a i
                if(proba[i] > MaxProb1){
                   MaxProb1 = proba[i]; 
                }
                
                // Auxiliar para o Log da Probabilidade Maxima
                if ( MaxProb1 != 0 )
                    Aux1 = (Aux3) * (Math.log(Pt) / Math.log (MaxProb1));
                else
                    Aux1 = 0;
                
                // Max da probabilidade de t+1 a 255
                MaxProb2 = 0;
                for (c=i+1; c<256; c++)
                {
                    if (proba[c] > MaxProb2 )	
                            MaxProb2 = proba[c];
                }
                
                // Auxiliar para o Log da Probabilidade Maxima
                if ( MaxProb2 != 0 )
                    Aux2 = (1 - Aux3) * (Math.log(1 - Pt) / Math.log (MaxProb2) );
                else
                    Aux2 = 0;

                // Calculo Final
                ni[i] = Aux1 + Aux2;
            }
        }
        
        max = 0;
	for (i = 0; i < 256; i++)
	{
            if (ni[i] > max)
            {
                max = ni[i];
                uiLimiar = i;
            }    
	}   
        
        // Visualizacao do valor de limiar de Otsu  em modo Debug
        // System.out.println(uiLimiar);
    
        //Cria a  imagem  binarizada 
	//Aloca a Matriz
	int [][] pBufferbinario = new int[Altura][Largura]; //Cria um PONTEIRO para a  imagem  binarizada 

        for( lin = 0; lin < Altura; lin++) {
            for( col = 0; col < Largura; col++) {
                Color x = new Color(img.getRGB(col,lin)); // leitura dos canais RGB de cada pixel
                // caso a imagem for colorida, gera uma imagem em ni­veis de cinza = (R+G+B)/3
		cinza = (int)((x.getGreen() + x.getRed() + x.getBlue())/3);
			
		if (cinza < uiLimiar)   pBufferbinario[lin][col] = 0;
		else    pBufferbinario[lin][col] = 1;  //multiplicado depois por 255
            }
	}

        //Aqui Gera a  imagem binaria 
        for(lin = 0; lin < Altura; lin++) {
		for(col = 0; col < Largura; col++) {
                    int atual = pBufferbinario[lin][col]* 255; //multiplicacao da img binaria por 255
                    Color novo = new Color(atual, atual, atual);
                    res.setRGB(col,lin, novo.getRGB()); // gravacao do valor binarizado para cada pixel
		}
        }

    return res;

}
    
    // Código de conversão  usando sRGB e não diretamente RGB
    private static double[] RGBtoLab(BufferedImage img, int image_x, int image_y){ 

        // exemplo de alocacao de vetor contendo as variáveis do espaço de cor Lab
        double[] pSaida = new double[3]; //Cria um PONTEIRO para armazenar cálculos
        
        Color orig = new Color(img.getRGB(image_x, image_y));  // lendo o valor RGB do pixel(i,j)
				
	// Conventendo RGB para sRGB - normalizando entre 0 e 1 
	// convert 0..255 into 0..1 
	double r = orig.getRed() / 255.0; 
        double g = orig.getGreen() /255.0;
	double b = orig.getBlue() / 255.0; 

	// Usando sRGB 
	if (r <= 0.04045)  r = r / 12.92; 
	else  r = Math.pow(((r + 0.055) / 1.055), 2.4); 
			
	if (g <= 0.04045) 	   g = g / 12.92; 
	else 	   g = Math.pow(((g + 0.055) / 1.055), 2.4); 

	if (b <= 0.04045) 	   b = b / 12.92; 
	else   b = Math.pow(((b + 0.055) / 1.055), 2.4); 

	r *= 100.0; // normalizando entre 0 e 100
        g *= 100.0; // normalizando entre 0 e 100
	b *= 100.0; // normalizando entre 0 e 100

        //Conversão de sRGB para Espaço de cor XYZ
        double Xx = r * 0.412424 + g * 0.357579 + b * 0.180464;
        double Yx = r * 0.212656 + g * 0.715158 + b * 0.0721856;
        double Zx = r * 0.0193324 + g * 0.119193 + b * 0.950444;

        //Conversão de XYZ para L* a* b* CORRIGIDO usando sRGB antes 
	// Usando Observer. = 2°, Illuminant = D65
	// D65 = {95.047, 100.000, 108.883};
	double x = Xx / 95.047; 
	double y = Yx / 100; 
	double z = Zx / 108.883; 

	if (x > 0.008856)    	x = Math.pow(x, 1.0/3.0);  
	else 	   				x = (7.787 * x) + (16.0/116.0); 
				 
	if (y > 0.008856) 	  	y = Math.pow(y, 1.0/3.0); 
	else 				  	y = (7.787 * y) + (16.0/116.0); 

	if (z > 0.008856) 	   z = Math.pow(z, 1.0/3.0); 
	else 				   z = (7.787 * z) + (16.0/116.0);  

	double L = (116.0 * y) - 16.0; 
	double a = 500.0 * (x - y); 
	double b_lab = 200.0 * (y - z); 
                
	//O retorno de cada indice do espaço de cores Lab para cada pixel da imagem dada
        pSaida[0] = L; 
        pSaida[1] = a; 
        pSaida[2] = b_lab; 
        
        return pSaida;
    }
    
    //
    public static BufferedImage EntropiaJohannsen(BufferedImage img){
        //Cria a  imagem  de sai­da  
        // Aloca a Matriz da imagem  de sai­da  
        BufferedImage res = new BufferedImage(img.getWidth(), img.getHeight(), img.getType());
        //Armazenando o tamanho da imagem
        int Largura = img.getWidth();    // Largura da imagem
        int Altura = img.getHeight();    // Altura da imagem
        int col, lin, i, cinza;
        double totalPixel = (double) Largura*Altura;
        double [] proba = new double[256];
        int [] histogram = new int[256];
        double [] sb = new double[256];
        double [] sw = new double[256];
        double [] tl = new double[256];
        int t;
        double sum1;
        double sum2;
        double minimo;
        
        // Inicializacao variaveis
   	int k, uiLimiar;        
        // inicializacao do Histograma
        for(i = 0; i < 256; i++) histogram[i]= 0;
        
        // Passo 1: calculo do Histograma  
	for( lin = 0; lin < Altura; lin++) {
            for( col = 0; col < Largura; col++) {
                Color x = new Color(img.getRGB(col,lin)); // leitura dos canais RGB de cada pixel
                // caso a imagem for colorida, gera uma imagem em ní­veis de cinza = (R+G+B)/3
                cinza = (int)((x.getGreen() + x.getRed() + x.getBlue())/3);		
                histogram[cinza]++; // Atualiza o histograma dos nÃ­veis de cinza (entre 0 e 255)
            }
	}
        
        //Passo 2: calculo das probabilidades a priori
        for (i = 0; i < 256; i++) {
            proba[i] = (double) ((histogram[i])/(double)(totalPixel));
            sb[i] = 0;
            sw[i] = 0;
            tl[i] = 0;
        }
        
        for (t = 0; t<256; t++) {
            sum1 = 0;
            sum2 = 0;
            
            for (int g=0; g<=t;g++) {
                sum1 = sum1 + proba[g];
            }
            for(int g=0;g<=t-1;g++) {
                sum2 = sum2 + proba[g];
            }
            
            if ((proba[t] != 0) && (sum1 != 0) && (sum2 != 0)) {
                sb[t] = (Math.log10((double)sum1)) - (1.0/sum1)*((proba[t]*Math.log10((double)proba[t]))) + (sum2*Math.log10((double)sum2));
            }
        }
        
        for (t = 0; t<256; t++) {
            sum1 = 0;
            sum2 = 0;
            
            for (int g = t; g < 256; g++) {
                sum1 = sum1 + proba[g];
            }
            for(int g = t+1; g < 256; g++) {
                sum2 = sum2 + proba[g];
            }
            
            if ((proba[t] != 0) && (sum1 != 0) && (sum2 != 0)) {
                sw[t] = (Math.log10((double)sum1)) - (1.0/sum1)*((proba[t]*Math.log10((double)proba[t]))) + (sum2*Math.log10((double)sum2));
            }           
        }
        // Calculo do Limiar de Johannsen
        for (int g = 0;g < 256; g++)
        {
            if ( (sb[g] != 0) && (sw[g] != 0) )
                tl[g] = sb[g] + sw[g];
        }

        int m=0;
        while(tl[m] == 0.0) {
            m++;
        }
        minimo = tl[m];
        uiLimiar = m;
        
        for(int g=m+1;g<256;g++) {
            if ( (tl[g] != 0) && (tl[g] < minimo) )
            {
                minimo = tl[g];
                 uiLimiar = g;
            }
	}

        // Visualizacao do valor de limiar de Otsu  em modo Debug
        System.out.println(uiLimiar);
        //Cria a  imagem  binarizada 
        //Aloca a Matriz
        int [][] pBufferbinario = new int[Altura][Largura]; //Cria um PONTEIRO para a  imagem  binarizada 

        for( lin = 0; lin < Altura; lin++) {
            for( col = 0; col < Largura; col++) {
                Color x = new Color(img.getRGB(col,lin)); // leitura dos canais RGB de cada pixel
                // caso a imagem for colorida, gera uma imagem em ni­veis de cinza = (R+G+B)/3
                cinza = (int)((x.getGreen() + x.getRed() + x.getBlue())/3);
                    
                if (cinza < uiLimiar)   pBufferbinario[lin][col] = 0;
                else pBufferbinario[lin][col] = 1;  //multiplicado depois por 255
            }
	}

        //Aqui Gera a  imagem binaria 
        for(lin = 0; lin < Altura; lin++) {
            for(col = 0; col < Largura; col++) {
                int atual = pBufferbinario[lin][col]* 255; //multiplicacao da img binaria por 255
                Color novo = new Color(atual, atual, atual);
                res.setRGB(col,lin, novo.getRGB()); // gravacao do valor binarizado para cada pixel
            }
        }

        return res;
        
    }
}