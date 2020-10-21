/*
 * Copyright (c) 2016. EAGER-CLI Alexander Peltzer
 * This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package Modules.preprocessing;

import IO.Communicator;
import Modules.AModule;
import com.google.common.io.Files;

import java.util.ArrayList;

/**
 * Created by peltzer on 23.01.14.
 */
public class ClipAndMerge extends AModule {

    public static final int DEFAULT = 0;
    private int currentConfiguration = DEFAULT;
    public static final int ADAPTER_CLIPPING_ONLY = 1;
    public static final int SINGLE_ENDED_ONLY = 2;



    public ClipAndMerge(Communicator c) {
        super(c);
    }

    public ClipAndMerge(Communicator c, int config){
        super(c);
        this.currentConfiguration = config;
    }

    @Override
    public void setParameters() {
        String output_stem = Files.getNameWithoutExtension(this.inputfile.get(0));
        this.outputfile = new ArrayList<String>();

        switch (currentConfiguration){
            case DEFAULT: this.parameters = getDefaultParameterList();
                appendMergedOnly();
                this.outputfile.add(getOutputfolder()+output_stem+".merged.fq.gz");
                break;
            case ADAPTER_CLIPPING_ONLY: this.parameters = getAdapterClippingOnlyParameterList();
                appendMergedOnly();
                this.outputfile.add(getOutputfolder()+output_stem+".clipped.fq.gz");
                break;
            case SINGLE_ENDED_ONLY: this.parameters = getSingleEndedOnlyParameterList();
                appendMergedOnly();

                this.outputfile.add(getOutputfolder()+output_stem+".fq.gz");
                break;


        }
       }


    private String[] getDefaultParameterList(){
        String output_stem = Files.getNameWithoutExtension(this.inputfile.get(0));
        return new String[]{"ClipAndMerge", "-in1",  getForwards(), "-in2", getReverse(),
                "-f", this.communicator.getMerge_fwadaptor(), "-r", this.communicator.getMerge_bwadaptor(),
                "-trim3p", this.communicator.getMerge_barcode3p(), "-trim5p", this.communicator.getMerge_barcode5p(),
                "-l", String.valueOf(this.communicator.getQuality_readlength()), "-m", String.valueOf(this.communicator.getMerge_min_adapter_overlap()), "-qt", "-q", String.valueOf(this.communicator.getQuality_minreadquality()),
                "-log", getOutputfolder()+ "stats.log",
                this.communicator.getMerge_advanced(),
                "-o ", getOutputfolder() +output_stem+".merged.fq.gz"};
    }

    private String[] getAdapterClippingOnlyParameterList(){
        String output_stem = Files.getNameWithoutExtension(this.inputfile.get(0));
        this.outputfile = new ArrayList<>();
        outputfile.add(getOutputfolder()+output_stem+".forwards.fq.gz");
        outputfile.add(getOutputfolder()+output_stem+".reverse.fq.gz");
        return new String[]{"ClipAndMerge", "-in1",  getForwards(), "-in2", getReverse(),
                "-f", this.communicator.getMerge_fwadaptor(), "-r", this.communicator.getMerge_bwadaptor(),
                "-l", String.valueOf(this.communicator.getQuality_readlength()),"-no_merging", "-qt", "-q", String.valueOf(this.communicator.getQuality_minreadquality()),
                "-log", getOutputfolder()+ "stats.log", "-m", String.valueOf(this.communicator.getMerge_min_adapter_overlap()),
                this.communicator.getMerge_advanced(),
                "-o ", getOutputfolder() +output_stem+".clipped.fq.gz", "-u", getOutputfolder()+output_stem+".forwards.fq.gz", getOutputfolder()+output_stem+".reverse.fq.gz"};

    }

    private String[] getSingleEndedOnlyParameterList() {
        String output_stem = Files.getNameWithoutExtension(this.inputfile.get(0));
        return new String[]{"ClipAndMerge", "-in1", getAll(),
                "-f", this.communicator.getMerge_fwadaptor(), "-r", this.communicator.getMerge_bwadaptor(),
                "-trim3p", this.communicator.getMerge_barcode3p(), "-trim5p", this.communicator.getMerge_barcode5p(),
                "-l", String.valueOf(this.communicator.getQuality_readlength()), "-qt", "-q", String.valueOf(this.communicator.getQuality_minreadquality()),
                "-log", getOutputfolder()+ "stats.log", "-m", String.valueOf(this.communicator.getMerge_min_adapter_overlap()),
                this.communicator.getMerge_advanced(),
                "-o ", getOutputfolder() +output_stem+".fq.gz"};
    }


    @Override
    public String getOutputfolder() {
        return  this.communicator.getGUI_resultspath()+"/1-AdapClip/";
    }

    @Override
    public String getModulename(){
        return super.getModulename() + getSubModuleName();
    }

    public String getSubModuleName(){
        switch(currentConfiguration){
            case DEFAULT: return "default";
            case ADAPTER_CLIPPING_ONLY: return "adapterclippingonly";
            case SINGLE_ENDED_ONLY: return "singleendedonly";
            default: return "default";
        }

    }




    private String getForwards(){

        String out = "";
        for(String s : this.communicator.getGUI_inputfiles()){
            if(s.contains("R1_") || s.contains("R1.")){
                out+= s + " ";
            }
        }
        return out;
    }

    private String getReverse(){
        String out = "";
        for(String s : this.communicator.getGUI_inputfiles()){
            if(s.contains("R2_") || s.contains("R2.")){
                out+= s + " ";
            }
        }
        return out;
    }

    private String getAll(){
        String out = "";
        for (String s : this.communicator.getGUI_inputfiles()){
            out+= s + " ";
        }
        return out;
    }

    private void appendMergedOnly(){
        String output_stem = Files.getNameWithoutExtension(this.inputfile.get(0));
        ArrayList<String> params = new ArrayList<String>();
        if(communicator.isMerge_keep_only_merged()) {
            for (String s : this.parameters) {
                params.add(s);
            }
            params.add("-u");
            params.add(getOutputfolder() + output_stem + ".forwards.unmerged.fq.gz");
            params.add(getOutputfolder() + output_stem + ".reverse.unmerged.fq.gz");


            this.parameters = new String[params.size()];

            for (int i = 0; i < params.size(); i++) {
                this.parameters[i] = params.get(i);
            }
        }

    }


}
