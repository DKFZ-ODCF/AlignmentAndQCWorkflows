package de.dkfz.b080.co;

import de.dkfz.roddy.plugins.BasePlugin;

/**

 * TODO Recreate class. Put in dependencies to other workflows, descriptions, capabilities (like ui settings, components) etc.
 */
public class QualityControlWorkflowPlugin extends BasePlugin {

<<<<<<< HEAD
    public static final String CURRENT_VERSION_STRING = "1.0.182";
    public static final String CURRENT_VERSION_BUILD_DATE = "Mon Dec 14 16:55:53 CET 2015";
=======
    public static final String CURRENT_VERSION_STRING = "1.0.177";
    public static final String CURRENT_VERSION_BUILD_DATE = "Fri Jul 01 16:46:08 CEST 2016";
>>>>>>> cbd96af... `bam` configuration value to provide an externally located BAM file as initial merged BAM into which to merge additional lane-BAMs. Bugfixes.

    @Override
    public String getVersionInfo() {
        return "Roddy plugin: " + this.getClass().getName() + ", V " + CURRENT_VERSION_STRING + " built at " + CURRENT_VERSION_BUILD_DATE;
    }
}
