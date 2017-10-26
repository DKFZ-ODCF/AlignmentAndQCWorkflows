package de.dkfz.b080.co;

import de.dkfz.roddy.plugins.BasePlugin;

/**

 * TODO Recreate class. Put in dependencies to other workflows, descriptions, capabilities (like ui settings, components) etc.
 */
public class QualityControlWorkflowPlugin extends BasePlugin {

	public static final String CURRENT_VERSION_STRING = "1.2.73";
	public static final String CURRENT_VERSION_BUILD_DATE = "Thu Oct 26 15:42:37 CEST 2017";

    @Override
    public String getVersionInfo() {
        return "Roddy plugin: " + this.getClass().getName() + ", V " + CURRENT_VERSION_STRING + " built at " + CURRENT_VERSION_BUILD_DATE;
    }
}
