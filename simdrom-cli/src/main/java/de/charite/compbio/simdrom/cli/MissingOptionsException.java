package de.charite.compbio.simdrom.cli;

import java.util.List;

import org.apache.commons.lang3.StringUtils;

public class MissingOptionsException extends Exception {

	/**
	 * 
	 */
	private static final long serialVersionUID = -2694951897378018172L;
	
	public MissingOptionsException(List<String> usedOptions, List<String> missingOptions) {
		super(buildMessage(usedOptions, missingOptions));

		
	}

	private static String buildMessage(List<String> usedOptions, List<String> missingOptions) {
		StringBuilder message = new StringBuilder("You cannot use the option(s) ");
		message.append(StringUtils.join(usedOptions, ", "));
		message.append(" without using option(s) ");
		message.append(StringUtils.join(usedOptions, ", "));
		return message.toString();
	}

}
