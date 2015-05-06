package de.charite.compbio.simdrom.cli;

import java.util.List;

import org.apache.commons.lang3.StringUtils;

public class NotAllowedCombinationOfOptionsException extends Exception {

	/**
	 * 
	 */
	private static final long serialVersionUID = -2694951897378018172L;
	
	public NotAllowedCombinationOfOptionsException(List<String> values) {
		super(buildMessage(values));

		
	}

	private static String buildMessage(List<String> values) {
		StringBuilder message = new StringBuilder("You cannot combine the options ");
		message.append(StringUtils.join(values, ", "));
		return message.toString();
	}

}
