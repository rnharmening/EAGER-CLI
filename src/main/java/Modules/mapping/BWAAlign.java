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

package Modules.mapping;

import IO.Communicator;
import Modules.AModule;
import com.google.common.io.Files;

/**
 * Created by peltzer on 24.01.14.
 */
public class BWAAlign extends AModule {
    public static final int DEFAULT = 0;
    private int currentconfiguration = DEFAULT;
    public static final int PAIREDENDWITHOUTMERGING = 1;
    public static final int MT = 2;
    public static final int MTREMAP = 3;

    public BWAAlign(Communicator c) {
        super(c);
    }

    public BWAAlign(Communicator c, int currentConfiguration) {
        super(c);
        this.currentconfiguration = currentConfiguration;
    }

    @Override
    public void setParameters() {
        switch (currentconfiguration) {
            case DEFAULT:
                this.parameters = getDefaultParameters();
                break;
            case PAIREDENDWITHOUTMERGING:
                this.parameters = getPairedEndWithoutMerging();
                break;
            case MT:
                this.parameters = getMTBWAAlignParameters();
                break;
        }
        this.outputfile = this.inputfile;
    }

    private String[] getDefaultParameters() {
        String output_stem = Files.getNameWithoutExtension(this.inputfile.get(0));
        String[] parameters = new String[]{"bwa", "aln", "-t",
                this.communicator.getCpucores(),
                this.communicator.getGUI_reference(), this.inputfile.get(0),
                "-n", this.communicator.getMapper_mismatches(),
                "-l", this.communicator.getMapper_seedlength(),
                this.communicator.getMapper_advanced(), "-f",
                getOutputfolder() + "/" + output_stem + ".sai"};

        return parameters;
    }

    private String[] getPairedEndWithoutMerging() {
        String output_stem0 = Files.getNameWithoutExtension(this.inputfile.get(0));
        String output_stem1 = Files.getNameWithoutExtension(this.inputfile.get(1));

        String commandOne = "bwa aln -t " + this.communicator.getCpucores() +
                " " + this.communicator.getGUI_reference() + " " + this.inputfile.get(0) + " " +
                "-n " + this.communicator.getMapper_mismatches() + " -l " + this.communicator.getMapper_seedlength() +
                " " + this.communicator.getMapper_advanced() + " -f " + getOutputfolder() + "/" + output_stem0 + ".sai";

        String commandTwo = "bwa aln -t " + this.communicator.getCpucores() +
                " " + this.communicator.getGUI_reference() + " " + this.inputfile.get(1) + " " +
                "-n " + this.communicator.getMapper_mismatches() + " -l " + this.communicator.getMapper_seedlength() +
                " " + this.communicator.getMapper_advanced() + " -f " + getOutputfolder() + "/" + output_stem1 + ".sai";

        String[] params = new String[]{"/bin/sh", "-c", commandOne + " && " + commandTwo};

        return params;
    }

    private String[] getMTBWAAlignParameters() {
        String output_stem = Files.getNameWithoutExtension(this.inputfile.get(0));
        return new String[]{"bwa", "aln", "-t",
                this.communicator.getCpucores(),
                this.communicator.getCM_referencemt_elong(), this.inputfile.get(0),
                "-n", this.communicator.getMapper_mismatches(),
                "-l", this.communicator.getMapper_seedlength(),
                this.communicator.getMapper_advanced(), "-f",
                getOutputfolder() + "/" + output_stem + ".MT.sai"};
    }


    @Override
    public String getOutputfolder() {
        switch (currentconfiguration) {
            case DEFAULT:
                return this.communicator.getGUI_resultspath() + "/3-Mapper";
            case MT:
                return this.communicator.getGUI_resultspath() + "/3-Mapper";
            default:
                return this.communicator.getGUI_resultspath() + "/3-Mapper";
        }
    }


}



