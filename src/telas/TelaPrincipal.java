/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package telas;

import java.awt.Image;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.imageio.ImageIO;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.filechooser.FileNameExtensionFilter;
import processamento.ProcessamentoImagem;

/**
 *
 * @author lucas
 */
public class TelaPrincipal extends javax.swing.JFrame {
    
    public File arqImagemOriginal;
    public BufferedImage imagemAtual;
    public boolean imagemFoiCarregada;

    /**
     * Creates new form TelaPrincipal
     */
    public TelaPrincipal() {
        initComponents();        
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        jMenu4 = new javax.swing.JMenu();
        jMenuItem5 = new javax.swing.JMenuItem();
        jMenuItem12 = new javax.swing.JMenuItem();
        jMenuItem14 = new javax.swing.JMenuItem();
        jMenuItem15 = new javax.swing.JMenuItem();
        jInternalFrame1 = new javax.swing.JInternalFrame();
        jScrollPane1 = new javax.swing.JScrollPane();
        jLabel1 = new javax.swing.JLabel();
        jMenuBar1 = new javax.swing.JMenuBar();
        jMenu1 = new javax.swing.JMenu();
        jMenuItem1 = new javax.swing.JMenuItem();
        jMenuItem2 = new javax.swing.JMenuItem();
        jMenuItem9 = new javax.swing.JMenuItem();
        jMenuItem3 = new javax.swing.JMenuItem();
        jMenu5 = new javax.swing.JMenu();
        jMenuItem18 = new javax.swing.JMenuItem();
        jMenuItem20 = new javax.swing.JMenuItem();
        jMenu2 = new javax.swing.JMenu();
        jMenu6 = new javax.swing.JMenu();
        jMenuItem6 = new javax.swing.JMenuItem();
        jMenuItem7 = new javax.swing.JMenuItem();
        jMenu3 = new javax.swing.JMenu();
        jMenuItem8 = new javax.swing.JMenuItem();
        jMenuItem10 = new javax.swing.JMenuItem();
        jMenu7 = new javax.swing.JMenu();
        jMenuItem17 = new javax.swing.JMenuItem();
        jMenuItem19 = new javax.swing.JMenuItem();
        jMenu11 = new javax.swing.JMenu();
        jMenuItem26 = new javax.swing.JMenuItem();
        jMenuItem27 = new javax.swing.JMenuItem();
        jMenu12 = new javax.swing.JMenu();
        jMenuItem28 = new javax.swing.JMenuItem();
        jMenuItem29 = new javax.swing.JMenuItem();
        jMenu13 = new javax.swing.JMenu();
        jMenuItem30 = new javax.swing.JMenuItem();
        jMenuItem31 = new javax.swing.JMenuItem();
        jMenu14 = new javax.swing.JMenu();
        jMenuItem32 = new javax.swing.JMenuItem();
        jMenuItem33 = new javax.swing.JMenuItem();
        jMenu9 = new javax.swing.JMenu();
        jMenuItem13 = new javax.swing.JMenuItem();
        jMenuItem4 = new javax.swing.JMenuItem();

        jMenu4.setText("jMenu4");

        jMenuItem5.setText("jMenuItem5");

        jMenuItem12.setText("jMenuItem12");

        jMenuItem14.setText("jMenuItem14");

        jMenuItem15.setText("jMenuItem15");

        jInternalFrame1.setVisible(true);

        javax.swing.GroupLayout jInternalFrame1Layout = new javax.swing.GroupLayout(jInternalFrame1.getContentPane());
        jInternalFrame1.getContentPane().setLayout(jInternalFrame1Layout);
        jInternalFrame1Layout.setHorizontalGroup(
            jInternalFrame1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 0, Short.MAX_VALUE)
        );
        jInternalFrame1Layout.setVerticalGroup(
            jInternalFrame1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 0, Short.MAX_VALUE)
        );

        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);
        setTitle("Detecção por Índice de Redness - Álef e Ketson");
        setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));
        setExtendedState(MAXIMIZED_BOTH);

        jLabel1.setBackground(new java.awt.Color(254, 254, 254));
        jLabel1.setVerticalAlignment(javax.swing.SwingConstants.TOP);
        jScrollPane1.setViewportView(jLabel1);

        jMenuBar1.setName(""); // NOI18N

        jMenu1.setText("Arquivo");

        jMenuItem1.setText("Carregar imagem");
        jMenuItem1.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItem1ActionPerformed(evt);
            }
        });
        jMenu1.add(jMenuItem1);

        jMenuItem2.setText("Restaurar original");
        jMenuItem2.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItem2ActionPerformed(evt);
            }
        });
        jMenu1.add(jMenuItem2);

        jMenuItem9.setText("Salvar");
        jMenuItem9.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItem9ActionPerformed(evt);
            }
        });
        jMenu1.add(jMenuItem9);

        jMenuItem3.setText("Sair");
        jMenuItem3.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItem3ActionPerformed(evt);
            }
        });
        jMenu1.add(jMenuItem3);

        jMenuBar1.add(jMenu1);

        jMenu5.setText("Índice de Redness");

        jMenuItem18.setText("Imagem carregada");
        jMenuItem18.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItem18ActionPerformed(evt);
            }
        });
        jMenu5.add(jMenuItem18);

        jMenuItem20.setText("Em lote");
        jMenuItem20.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItem20ActionPerformed(evt);
            }
        });
        jMenu5.add(jMenuItem20);

        jMenuBar1.add(jMenu5);

        jMenu2.setText("Técnicas de Binarização");

        jMenu6.setText("Binarização de Otsu");

        jMenuItem6.setText("Imagem carregada");
        jMenuItem6.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItem6ActionPerformed(evt);
            }
        });
        jMenu6.add(jMenuItem6);

        jMenuItem7.setText("Em lote");
        jMenuItem7.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItem7ActionPerformed(evt);
            }
        });
        jMenu6.add(jMenuItem7);

        jMenu2.add(jMenu6);

        jMenu3.setText("Entropia de Pun");

        jMenuItem8.setText("Imagem carregada");
        jMenuItem8.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItem8ActionPerformed(evt);
            }
        });
        jMenu3.add(jMenuItem8);

        jMenuItem10.setText("Em lote");
        jMenuItem10.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItem10ActionPerformed(evt);
            }
        });
        jMenu3.add(jMenuItem10);

        jMenu2.add(jMenu3);

        jMenu7.setText("Entropia Johannsen");

        jMenuItem17.setText("Imagem carregada");
        jMenuItem17.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItem17ActionPerformed(evt);
            }
        });
        jMenu7.add(jMenuItem17);

        jMenuItem19.setText("Em lote");
        jMenuItem19.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItem19ActionPerformed(evt);
            }
        });
        jMenu7.add(jMenuItem19);

        jMenu2.add(jMenu7);

        jMenuBar1.add(jMenu2);

        jMenu11.setText("Índice de Redness - Algoritmo 1");

        jMenuItem26.setText("Imagem carregada");
        jMenuItem26.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItem26ActionPerformed(evt);
            }
        });
        jMenu11.add(jMenuItem26);

        jMenuItem27.setText("Em lote");
        jMenuItem27.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItem27ActionPerformed(evt);
            }
        });
        jMenu11.add(jMenuItem27);

        jMenuBar1.add(jMenu11);

        jMenu12.setText("Índice de Redness - Algoritmo 2");

        jMenuItem28.setText("Imagem carregada");
        jMenuItem28.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItem28ActionPerformed(evt);
            }
        });
        jMenu12.add(jMenuItem28);

        jMenuItem29.setText("Em lote");
        jMenuItem29.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItem29ActionPerformed(evt);
            }
        });
        jMenu12.add(jMenuItem29);

        jMenuBar1.add(jMenu12);

        jMenu13.setText("Índice de Redness - Fusão");

        jMenuItem30.setText("Imagem carregada");
        jMenuItem30.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItem30ActionPerformed(evt);
            }
        });
        jMenu13.add(jMenuItem30);

        jMenuItem31.setText("Em lote");
        jMenuItem31.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItem31ActionPerformed(evt);
            }
        });
        jMenu13.add(jMenuItem31);

        jMenuBar1.add(jMenu13);

        jMenu14.setText("DETECÇÃO FINAL");

        jMenuItem32.setText("Imagem carregada");
        jMenuItem32.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItem32ActionPerformed(evt);
            }
        });
        jMenu14.add(jMenuItem32);

        jMenuItem33.setText("Em lote");
        jMenuItem33.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItem33ActionPerformed(evt);
            }
        });
        jMenu14.add(jMenuItem33);

        jMenuBar1.add(jMenu14);

        jMenu9.setText("Sobre");

        jMenuItem13.setText("Autor da Interface");
        jMenuItem13.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItem13ActionPerformed(evt);
            }
        });
        jMenu9.add(jMenuItem13);

        jMenuItem4.setText("Autor do Programa de Detecção");
        jMenuItem4.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItem4ActionPerformed(evt);
            }
        });
        jMenu9.add(jMenuItem4);

        jMenuBar1.add(jMenu9);

        setJMenuBar(jMenuBar1);

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jScrollPane1, javax.swing.GroupLayout.DEFAULT_SIZE, 1225, Short.MAX_VALUE)
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addComponent(jScrollPane1, javax.swing.GroupLayout.DEFAULT_SIZE, 611, Short.MAX_VALUE)
                .addContainerGap())
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void jMenuItem1ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItem1ActionPerformed
        JFileChooser fs = new JFileChooser(new File("/home/"));
        fs.setFileFilter(new FileNameExtensionFilter("Arquivos de Imagem (png, jpg, jpeg, bmp)","png","jpg","jpeg","bmp"));
        fs.setAcceptAllFileFilterUsed(false);
        fs.setDialogTitle("Selecione o arquivo de imagem...");
        int retorno = fs.showOpenDialog(this);        
        if(retorno == JFileChooser.APPROVE_OPTION){
            try {
                arqImagemOriginal = fs.getSelectedFile();
                imagemAtual = ImageIO.read(arqImagemOriginal);
                jLabel1.setIcon(new javax.swing.ImageIcon(imagemAtual));
                imagemFoiCarregada = true;
            } catch (IOException ex) {
                Logger.getLogger(TelaPrincipal.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }//GEN-LAST:event_jMenuItem1ActionPerformed

    private void jMenuItem9ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItem9ActionPerformed
        JFrame f=new JFrame();
        if(imagemFoiCarregada){
            JFileChooser fs = new JFileChooser(new File(arqImagemOriginal.getParent()));
            fs.setFileFilter(new FileNameExtensionFilter("PNG Images","png"));
            fs.setDialogTitle("Escolha onde salvar a imagem posterizada:");
            int retorno = fs.showSaveDialog(this);
            if(retorno ==JFileChooser.APPROVE_OPTION){   
                 File fo = new File(fs.getSelectedFile().getPath()+"_RedCIRG2"+".png");
                try {
                    ImageIO.write(imagemAtual, "png", fo);
                    JOptionPane.showMessageDialog(f, "Imagem salva com sucesso", "Mensagem", JOptionPane.INFORMATION_MESSAGE);
                } catch (IOException ex) {
                    JOptionPane.showMessageDialog(f, "Não foi possível salvar a imagem!", "Aviso", JOptionPane.WARNING_MESSAGE);
                    Logger.getLogger(TelaPrincipal.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
        }else{
            JOptionPane.showMessageDialog(f, "A imagem não foi carregada!", "Aviso", JOptionPane.WARNING_MESSAGE);
        }
    }//GEN-LAST:event_jMenuItem9ActionPerformed

    private void jMenuItem2ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItem2ActionPerformed
        if(imagemFoiCarregada){    
            try {
                imagemAtual = ImageIO.read(arqImagemOriginal);
            } catch (IOException ex) {
                Logger.getLogger(TelaPrincipal.class.getName()).log(Level.SEVERE, null, ex);
            }
            jLabel1.setIcon(new javax.swing.ImageIcon(imagemAtual));
        }
    }//GEN-LAST:event_jMenuItem2ActionPerformed

    private void jMenuItem13ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItem13ActionPerformed
        JFrame f=new JFrame();
        JOptionPane.showMessageDialog(f, "Interface desenvolvida por Álef e Ketson para a disciplina de PDI 2023/1, tendo sido ministrada pelo Prof. Jacques Facon.", "Autor", JOptionPane.INFORMATION_MESSAGE);
    }//GEN-LAST:event_jMenuItem13ActionPerformed
    
    private void jMenuItem3ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItem3ActionPerformed
        JOptionPane.showMessageDialog(null, "Volte sempre!", "Mensagem", JOptionPane.INFORMATION_MESSAGE);
        System.exit(0);
    }//GEN-LAST:event_jMenuItem3ActionPerformed

    private void jMenuItem18ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItem18ActionPerformed
        JFrame f=new JFrame();
        if(imagemFoiCarregada){
            try{
                BufferedImage imagemPosterizada = ProcessamentoImagem.detectarPeleHumana(imagemAtual); 
                if(imagemPosterizada == null) {
                    JOptionPane.showMessageDialog(f, "Não foi possível realizar a posterizacão!", "Aviso", JOptionPane.WARNING_MESSAGE);
                }else{
                    imagemAtual = imagemPosterizada;
                    jLabel1.setIcon(new javax.swing.ImageIcon(imagemAtual));
                }
            }catch( NullPointerException | IllegalArgumentException ex){
                JOptionPane.showMessageDialog(f, "Não foi possível realizar a posterizacão!", "Aviso", JOptionPane.WARNING_MESSAGE);
                Logger.getLogger(TelaPrincipal.class.getName()).log(Level.SEVERE, null, ex);
            }
        }else{
            JOptionPane.showMessageDialog(f, "A imagem não foi carregada!", "Aviso", JOptionPane.WARNING_MESSAGE);
        }
    }//GEN-LAST:event_jMenuItem18ActionPerformed

    private void jMenuItem20ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItem20ActionPerformed
        TelaDetectaPelePasta TDPP = new TelaDetectaPelePasta();
        TDPP.setLocationRelativeTo(null);
        TDPP.setVisible(true);
    }//GEN-LAST:event_jMenuItem20ActionPerformed

    private void jMenuItem4ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItem4ActionPerformed
        JFrame f=new JFrame();
        JOptionPane.showMessageDialog(f, "Programa desenvolvido por Álef e Ketson para a disciplina de PDI 2022/2, tendo sido ministrada pelo Prof. Jacques Facon.", "Autor", JOptionPane.INFORMATION_MESSAGE);
    }//GEN-LAST:event_jMenuItem4ActionPerformed

    private void jMenuItem6ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItem6ActionPerformed
        JFrame f=new JFrame();
        if(imagemFoiCarregada){
            try{
                BufferedImage imagemPosterizada = ProcessamentoImagem.BinarizacaoOtsu(imagemAtual); 
                if(imagemPosterizada == null) {
                    JOptionPane.showMessageDialog(f, "Não foi possível realizar a posterizacão!", "Aviso", JOptionPane.WARNING_MESSAGE);
                }else{
                    imagemAtual = imagemPosterizada;
                    jLabel1.setIcon(new javax.swing.ImageIcon(imagemAtual));
                }
            }catch( NullPointerException | IllegalArgumentException ex){
                JOptionPane.showMessageDialog(f, "Não foi possível realizar a posterizacão!", "Aviso", JOptionPane.WARNING_MESSAGE);
                Logger.getLogger(TelaPrincipal.class.getName()).log(Level.SEVERE, null, ex);
            }
        }else{
            JOptionPane.showMessageDialog(f, "A imagem não foi carregada!", "Aviso", JOptionPane.WARNING_MESSAGE);
        }
    }//GEN-LAST:event_jMenuItem6ActionPerformed

    private void jMenuItem7ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItem7ActionPerformed
        TelaBinarizacaoOtsu TBO = new TelaBinarizacaoOtsu();
        TBO.setLocationRelativeTo(null);
        TBO.setVisible(true);
    }//GEN-LAST:event_jMenuItem7ActionPerformed

    private void jMenuItem8ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItem8ActionPerformed
        JFrame f=new JFrame();
        if(imagemFoiCarregada){
            try{
                BufferedImage imagemPosterizada = ProcessamentoImagem.EntropiaPun(imagemAtual); 
                if(imagemPosterizada == null) {
                    JOptionPane.showMessageDialog(f, "Não foi possível realizar a posterizacão!", "Aviso", JOptionPane.WARNING_MESSAGE);
                }else{
                    imagemAtual = imagemPosterizada;
                    jLabel1.setIcon(new javax.swing.ImageIcon(imagemAtual));
                }
            }catch( NullPointerException | IllegalArgumentException ex){
                JOptionPane.showMessageDialog(f, "Não foi possível realizar a posterizacão!", "Aviso", JOptionPane.WARNING_MESSAGE);
                Logger.getLogger(TelaPrincipal.class.getName()).log(Level.SEVERE, null, ex);
            }
        }else{
            JOptionPane.showMessageDialog(f, "A imagem não foi carregada!", "Aviso", JOptionPane.WARNING_MESSAGE);
        }
    }//GEN-LAST:event_jMenuItem8ActionPerformed

    private void jMenuItem10ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItem10ActionPerformed
        TelaEntropiaPun TBO = new TelaEntropiaPun();
        TBO.setLocationRelativeTo(null);
        TBO.setVisible(true);
    }//GEN-LAST:event_jMenuItem10ActionPerformed

    private void jMenuItem26ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItem26ActionPerformed
        JFrame f=new JFrame();
        if(imagemFoiCarregada){
            try{
                BufferedImage imagemPosterizada = ProcessamentoImagem.rednessDetection1(imagemAtual);
                if(imagemPosterizada == null) {
                    JOptionPane.showMessageDialog(f, "Não foi possível realizar a posterizacão!", "Aviso", JOptionPane.WARNING_MESSAGE);
                }else{
                    imagemAtual = imagemPosterizada;
                    jLabel1.setIcon(new javax.swing.ImageIcon(imagemAtual));
                }
            }catch( NullPointerException | IllegalArgumentException ex){
                JOptionPane.showMessageDialog(f, "Não foi possível realizar a posterizacão!", "Aviso", JOptionPane.WARNING_MESSAGE);
                Logger.getLogger(TelaPrincipal.class.getName()).log(Level.SEVERE, null, ex);
            }
        }else{
            JOptionPane.showMessageDialog(f, "A imagem não foi carregada!", "Aviso", JOptionPane.WARNING_MESSAGE);
        }
    }//GEN-LAST:event_jMenuItem26ActionPerformed

    private void jMenuItem27ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItem27ActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_jMenuItem27ActionPerformed

    private void jMenuItem28ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItem28ActionPerformed
        JFrame f=new JFrame();
        if(imagemFoiCarregada){
            try{
                BufferedImage imagemPosterizada = ProcessamentoImagem.rednessDetection2(imagemAtual);
                if(imagemPosterizada == null) {
                    JOptionPane.showMessageDialog(f, "Não foi possível realizar a posterizacão!", "Aviso", JOptionPane.WARNING_MESSAGE);
                }else{
                    imagemAtual = imagemPosterizada;
                    jLabel1.setIcon(new javax.swing.ImageIcon(imagemAtual));
                }
            }catch( NullPointerException | IllegalArgumentException ex){
                JOptionPane.showMessageDialog(f, "Não foi possível realizar a posterizacão!", "Aviso", JOptionPane.WARNING_MESSAGE);
                Logger.getLogger(TelaPrincipal.class.getName()).log(Level.SEVERE, null, ex);
            }
        }else{
            JOptionPane.showMessageDialog(f, "A imagem não foi carregada!", "Aviso", JOptionPane.WARNING_MESSAGE);
        }
    }//GEN-LAST:event_jMenuItem28ActionPerformed

    private void jMenuItem29ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItem29ActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_jMenuItem29ActionPerformed

    private void jMenuItem30ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItem30ActionPerformed
        JFrame f=new JFrame();
        if(imagemFoiCarregada){
            try{
                BufferedImage imagemPosterizada = ProcessamentoImagem.rednessFusion(imagemAtual);
                if(imagemPosterizada == null) {
                    JOptionPane.showMessageDialog(f, "Não foi possível realizar a posterizacão!", "Aviso", JOptionPane.WARNING_MESSAGE);
                }else{
                    imagemAtual = imagemPosterizada;
                    jLabel1.setIcon(new javax.swing.ImageIcon(imagemAtual));
                }
            }catch( NullPointerException | IllegalArgumentException ex){
                JOptionPane.showMessageDialog(f, "Não foi possível realizar a posterizacão!", "Aviso", JOptionPane.WARNING_MESSAGE);
                Logger.getLogger(TelaPrincipal.class.getName()).log(Level.SEVERE, null, ex);
            }
        }else{
            JOptionPane.showMessageDialog(f, "A imagem não foi carregada!", "Aviso", JOptionPane.WARNING_MESSAGE);
        }
    }//GEN-LAST:event_jMenuItem30ActionPerformed

    private void jMenuItem31ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItem31ActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_jMenuItem31ActionPerformed

    private void jMenuItem32ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItem32ActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_jMenuItem32ActionPerformed

    private void jMenuItem33ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItem33ActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_jMenuItem33ActionPerformed

    private void jMenuItem17ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItem17ActionPerformed
        JFrame f=new JFrame();
        if(imagemFoiCarregada){
            try{
                BufferedImage imagemPosterizada = ProcessamentoImagem.EntropiaJohannsen(imagemAtual); 
                if(imagemPosterizada == null) {
                    JOptionPane.showMessageDialog(f, "Não foi possível realizar a posterizacão!", "Aviso", JOptionPane.WARNING_MESSAGE);
                }else{
                    imagemAtual = imagemPosterizada;
                    jLabel1.setIcon(new javax.swing.ImageIcon(imagemAtual));
                }
            }catch( NullPointerException | IllegalArgumentException ex){
                JOptionPane.showMessageDialog(f, "Não foi possível realizar a posterizacão!", "Aviso", JOptionPane.WARNING_MESSAGE);
                Logger.getLogger(TelaPrincipal.class.getName()).log(Level.SEVERE, null, ex);
            }
        }else{
            JOptionPane.showMessageDialog(f, "A imagem não foi carregada!", "Aviso", JOptionPane.WARNING_MESSAGE);
        }
    }//GEN-LAST:event_jMenuItem17ActionPerformed

    private void jMenuItem19ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jMenuItem19ActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_jMenuItem19ActionPerformed

    /**
     * @param args the command line arguments
     */
    public static void main(String args[]) {
        /* Set the Nimbus look and feel */
        //<editor-fold defaultstate="collapsed" desc=" Look and feel setting code (optional) ">
        /* If Nimbus (introduced in Java SE 6) is not available, stay with the default look and feel.
         * For details see http://download.oracle.com/javase/tutorial/uiswing/lookandfeel/plaf.html 
         */
        try {
            for (javax.swing.UIManager.LookAndFeelInfo info : javax.swing.UIManager.getInstalledLookAndFeels()) {
                if ("FlatLaf".equals(info.getName())) {
                    javax.swing.UIManager.setLookAndFeel(info.getClassName());
                    break;
                }
            }
        } catch (ClassNotFoundException ex) {
            java.util.logging.Logger.getLogger(TelaPrincipal.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (InstantiationException ex) {
            java.util.logging.Logger.getLogger(TelaPrincipal.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (IllegalAccessException ex) {
            java.util.logging.Logger.getLogger(TelaPrincipal.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (javax.swing.UnsupportedLookAndFeelException ex) {
            java.util.logging.Logger.getLogger(TelaPrincipal.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        }
        //</editor-fold>

        /* Create and display the form */
        java.awt.EventQueue.invokeLater(new Runnable() {
            public void run() {
                new TelaPrincipal().setVisible(true);
            }
        });
        
        
    }

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JInternalFrame jInternalFrame1;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JMenu jMenu1;
    private javax.swing.JMenu jMenu11;
    private javax.swing.JMenu jMenu12;
    private javax.swing.JMenu jMenu13;
    private javax.swing.JMenu jMenu14;
    private javax.swing.JMenu jMenu2;
    private javax.swing.JMenu jMenu3;
    private javax.swing.JMenu jMenu4;
    private javax.swing.JMenu jMenu5;
    private javax.swing.JMenu jMenu6;
    private javax.swing.JMenu jMenu7;
    private javax.swing.JMenu jMenu9;
    private javax.swing.JMenuBar jMenuBar1;
    private javax.swing.JMenuItem jMenuItem1;
    private javax.swing.JMenuItem jMenuItem10;
    private javax.swing.JMenuItem jMenuItem12;
    private javax.swing.JMenuItem jMenuItem13;
    private javax.swing.JMenuItem jMenuItem14;
    private javax.swing.JMenuItem jMenuItem15;
    private javax.swing.JMenuItem jMenuItem17;
    private javax.swing.JMenuItem jMenuItem18;
    private javax.swing.JMenuItem jMenuItem19;
    private javax.swing.JMenuItem jMenuItem2;
    private javax.swing.JMenuItem jMenuItem20;
    private javax.swing.JMenuItem jMenuItem26;
    private javax.swing.JMenuItem jMenuItem27;
    private javax.swing.JMenuItem jMenuItem28;
    private javax.swing.JMenuItem jMenuItem29;
    private javax.swing.JMenuItem jMenuItem3;
    private javax.swing.JMenuItem jMenuItem30;
    private javax.swing.JMenuItem jMenuItem31;
    private javax.swing.JMenuItem jMenuItem32;
    private javax.swing.JMenuItem jMenuItem33;
    private javax.swing.JMenuItem jMenuItem4;
    private javax.swing.JMenuItem jMenuItem5;
    private javax.swing.JMenuItem jMenuItem6;
    private javax.swing.JMenuItem jMenuItem7;
    private javax.swing.JMenuItem jMenuItem8;
    private javax.swing.JMenuItem jMenuItem9;
    private javax.swing.JScrollPane jScrollPane1;
    // End of variables declaration//GEN-END:variables
}
