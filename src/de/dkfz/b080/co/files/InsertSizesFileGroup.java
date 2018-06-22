/*
 * Copyright (c) 2018 German Cancer Research Center (DKFZ).
 *
 * Distributed under the MIT License (license terms are at https://github.com/TheRoddyWMS/AlignmentAndQCWorkflows).
 */

package de.dkfz.b080.co.files;

import de.dkfz.roddy.knowledge.files.BaseFile;
import de.dkfz.roddy.knowledge.files.FileGroup;

import java.util.List;

/**
 * @author michael
 */
public class InsertSizesFileGroup extends FileGroup {

    private final InsertSizesPlotFile plotFile;
    private final InsertSizesTextFile textFile;
    private final InsertSizesValueFile valueFile;

    public InsertSizesFileGroup(List<BaseFile> allFiles) {
        super(allFiles);
        InsertSizesTextFile tf = null;
        InsertSizesPlotFile pf = null;
        InsertSizesValueFile vf = null;
        for (BaseFile bf : allFiles) {
            if (bf instanceof InsertSizesTextFile) {
                tf = (InsertSizesTextFile) bf;
            } else if (bf instanceof InsertSizesPlotFile) {
                pf = (InsertSizesPlotFile) bf;
            } else if (bf instanceof InsertSizesValueFile) {
                vf = (InsertSizesValueFile) bf;
            }
        }
        this.plotFile = pf;
        this.textFile = tf;
        this.valueFile = vf;
    }

    public InsertSizesPlotFile getPlotFile() {
        return plotFile;
    }

    public InsertSizesTextFile getTextFile() {
        return textFile;
    }

    public InsertSizesValueFile getValueFile() {
        return valueFile;
    }
}
